#!/usr/bin/env python3
"""Create FASTA and sequence logo from CSV sequences aligned to a PDB reference chain."""

from __future__ import annotations

import argparse
import re
import sys
from collections import OrderedDict
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import seq1

try:
    import logomaker
except ImportError:  # pragma: no cover - runtime dependency check
    logomaker = None


AA_ALPHABET = list("ACDEFGHIKLMNPQRSTVWY")
VALID_AA_PATTERN = re.compile(r"^[A-Z]+$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Read sequences from a CSV column, add a reference sequence from a PDB chain, "
            "run MSA, write FASTA, and render a sequence logo aligned to the reference."
        )
    )
    parser.add_argument("--csv", required=True, help="Path to CSV file containing sequences")
    parser.add_argument("--pdb", required=True, help="Path to reference PDB file")
    parser.add_argument("--chain", required=True, help="Reference chain ID (e.g. A)")
    parser.add_argument(
        "--sequence-column",
        default="sequence",
        help="CSV column with amino-acid sequences (default: sequence)",
    )
    parser.add_argument(
        "--id-column",
        default="description",
        help="CSV column for FASTA IDs (default: description)",
    )
    parser.add_argument(
        "--outdir",
        default=".",
        help="Output directory for FASTA/logo files (default: current directory)",
    )
    parser.add_argument(
        "--prefix",
        default="csv_logo",
        help="Output filename prefix (default: csv_logo)",
    )
    parser.add_argument(
        "--tick-step",
        type=int,
        default=5,
        help="Tick spacing for alignment/reference x-axis labels (default: 5)",
    )
    parser.add_argument(
        "--title",
        default=None,
        help="Optional plot title (default: derived from input filenames)",
    )
    parser.add_argument(
        "--allow-gaps",
        action="store_true",
        help=(
            "Allow gapped alignment via built-in star alignment when sequences differ in length. "
            "Default behavior is strict equal-length, gap-free comparison."
        ),
    )
    return parser.parse_args()


def sanitize_fasta_id(raw: str, fallback: str) -> str:
    value = (raw or "").strip()
    if not value:
        value = fallback
    value = re.sub(r"\s+", "_", value)
    value = re.sub(r"[^A-Za-z0-9_.|:-]", "_", value)
    return value


def load_csv_sequences(
    csv_path: Path,
    sequence_column: str,
    id_column: str,
) -> tuple[pd.DataFrame, OrderedDict[str, str]]:
    if not csv_path.is_file():
        raise FileNotFoundError(f"CSV file not found: {csv_path}")

    df = pd.read_csv(csv_path)
    if sequence_column not in df.columns:
        raise ValueError(
            f"Sequence column '{sequence_column}' not found in CSV. "
            f"Available columns: {', '.join(df.columns)}"
        )

    if id_column not in df.columns:
        print(
            f"Warning: id column '{id_column}' not found, using row-index-based IDs.",
            file=sys.stderr,
        )

    seq_df = df.copy()
    seq_df[sequence_column] = (
        seq_df[sequence_column]
        .astype(str)
        .str.replace(" ", "", regex=False)
        .str.strip()
        .str.upper()
    )

    seq_df = seq_df[seq_df[sequence_column].notna()]
    seq_df = seq_df[seq_df[sequence_column] != ""]

    if seq_df.empty:
        raise ValueError("No non-empty sequences found in CSV after cleaning.")

    bad_rows = []
    for idx, seq in seq_df[sequence_column].items():
        if not VALID_AA_PATTERN.match(seq):
            bad_rows.append((idx, seq))

    if bad_rows:
        preview = "; ".join(f"row {idx}: {seq}" for idx, seq in bad_rows[:5])
        raise ValueError(
            "Found sequences with non-alphabetic characters. "
            f"Examples: {preview}"
        )

    records: OrderedDict[str, str] = OrderedDict()
    seen: dict[str, int] = {}

    for i, row in seq_df.iterrows():
        raw_id = str(row[id_column]) if id_column in seq_df.columns else ""
        base_id = sanitize_fasta_id(raw_id, fallback=f"seq_{i}")
        suffix = seen.get(base_id, 0)
        seen[base_id] = suffix + 1
        final_id = base_id if suffix == 0 else f"{base_id}_{suffix + 1}"
        records[final_id] = row[sequence_column]

    if not records:
        raise ValueError("No valid sequence records could be created from CSV.")

    return seq_df[[sequence_column]].copy(), records


def extract_reference_sequence(
    pdb_path: Path,
    chain_id: str,
) -> tuple[str, list[str]]:
    if not pdb_path.is_file():
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_path.stem, str(pdb_path))

    if 0 not in structure:
        raise ValueError(f"No model 0 found in PDB: {pdb_path}")

    model = structure[0]
    if chain_id not in [c.id for c in model]:
        raise ValueError(f"Chain '{chain_id}' not found in PDB: {pdb_path.name}")

    chain = model[chain_id]

    sequence_chars: list[str] = []
    residue_labels: list[str] = []

    for residue in chain:
        if not is_aa(residue, standard=True):
            continue
        aa = seq1(residue.resname, undef_code="X")
        if aa == "X":
            continue

        _, resseq, icode = residue.id
        icode = (icode or "").strip()
        label = f"{resseq}{icode}" if icode else str(resseq)

        sequence_chars.append(aa)
        residue_labels.append(label)

    if not sequence_chars:
        raise ValueError(
            f"No standard amino-acid residues found in chain '{chain_id}' of {pdb_path.name}"
        )

    return "".join(sequence_chars), residue_labels


def _parse_clustal_gapped(clustal_str: str) -> tuple[str, str]:
    parts = [[], []]
    for line in clustal_str.splitlines():
        for i, prefix in enumerate(("sequence_0", "sequence_1")):
            if line.startswith(prefix):
                parts[i].append(line.split()[-1])
    return "".join(parts[0]), "".join(parts[1])


def _build_aligner() -> PairwiseAligner:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -3
    aligner.extend_gap_score = -0.5
    return aligner


def star_align(sequences: OrderedDict[str, str]) -> OrderedDict[str, str]:
    names = list(sequences.keys())
    seqs = list(sequences.values())
    n = len(seqs)

    if n == 1:
        return sequences

    aligner = _build_aligner()

    scores = [sum(aligner.score(seqs[i], seqs[j]) for j in range(n) if j != i) for i in range(n)]
    center_idx = scores.index(max(scores))
    center = seqs[center_idx]
    center_len = len(center)

    ins_schedule = [{} for _ in range(n)]
    match_schedule = [{} for _ in range(n)]

    for i, seq in enumerate(seqs):
        aln = aligner.align(center, seq)[0]
        t_gapped, q_gapped = _parse_clustal_gapped(aln.format("clustal"))

        center_pos = 0
        pending_ins: list[str] = []
        for t_ch, q_ch in zip(t_gapped, q_gapped):
            if t_ch == "-":
                pending_ins.append(q_ch)
            else:
                if pending_ins:
                    existing = ins_schedule[i].get(center_pos, [])
                    ins_schedule[i][center_pos] = existing + pending_ins
                    pending_ins = []
                match_schedule[i][center_pos] = q_ch
                center_pos += 1
        if pending_ins:
            existing = ins_schedule[i].get(center_len, [])
            ins_schedule[i][center_len] = existing + pending_ins

    max_ins: dict[int, int] = {}
    for i in range(n):
        for pos, ins_list in ins_schedule[i].items():
            max_ins[pos] = max(max_ins.get(pos, 0), len(ins_list))

    final_seqs = ["" for _ in range(n)]
    for pos in range(center_len + 1):
        ins_width = max_ins.get(pos, 0)
        if ins_width > 0:
            for i in range(n):
                ins_chars = ins_schedule[i].get(pos, [])
                padded = "".join(ins_chars) + "-" * (ins_width - len(ins_chars))
                final_seqs[i] += padded
        if pos < center_len:
            for i in range(n):
                final_seqs[i] += match_schedule[i].get(pos, "-")

    return OrderedDict(zip(names, final_seqs))


def write_fasta(records: OrderedDict[str, str], output_path: Path) -> None:
    seq_records = [SeqRecord(Seq(seq), id=name, description="") for name, seq in records.items()]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(seq_records, str(output_path), "fasta")


def assert_equal_lengths(reference_sequence: str, design_records: OrderedDict[str, str]) -> None:
    if not design_records:
        raise ValueError("No design sequences were found in the CSV.")

    ref_len = len(reference_sequence)
    mismatched: list[tuple[str, int]] = []

    for name, seq in design_records.items():
        if len(seq) != ref_len:
            mismatched.append((name, len(seq)))

    if mismatched:
        preview = "; ".join(f"{name}={seq_len}" for name, seq_len in mismatched[:10])
        raise ValueError(
            "Reference and design sequences must be the same length for gap-free comparison. "
            f"Reference length={ref_len}; mismatches: {preview}"
        )


def compute_logo_matrix(aligned_sequences: list[str]) -> pd.DataFrame:
    if not aligned_sequences:
        raise ValueError("No aligned design sequences provided for logo computation.")

    aln_len = len(aligned_sequences[0])
    if any(len(seq) != aln_len for seq in aligned_sequences):
        raise ValueError("Aligned sequences have inconsistent lengths.")

    rows: list[dict[str, float]] = []
    alphabet_set = set(AA_ALPHABET)

    for col_idx in range(aln_len):
        counts = {aa: 0 for aa in AA_ALPHABET}
        observed = 0
        for seq in aligned_sequences:
            aa = seq[col_idx]
            if aa in alphabet_set:
                counts[aa] += 1
                observed += 1
        if observed > 0:
            rows.append({aa: counts[aa] / observed for aa in AA_ALPHABET})
        else:
            rows.append({aa: 0.0 for aa in AA_ALPHABET})

    return pd.DataFrame(rows)


def map_alignment_to_reference_labels(aligned_reference: str, ref_labels: list[str]) -> list[str | None]:
    mapped: list[str | None] = []
    ref_idx = 0

    for aa in aligned_reference:
        if aa == "-":
            mapped.append(None)
            continue

        if ref_idx >= len(ref_labels):
            mapped.append(None)
            continue

        mapped.append(ref_labels[ref_idx])
        ref_idx += 1

    return mapped


def _clustal_color_scheme() -> dict[str, str]:
    return {
        "A": "#f0a100",
        "C": "#f0a100",
        "F": "#f0a100",
        "I": "#f0a100",
        "L": "#f0a100",
        "M": "#f0a100",
        "V": "#f0a100",
        "W": "#f0a100",
        "Y": "#f0a100",
        "H": "#4040f2",
        "K": "#4040f2",
        "R": "#4040f2",
        "D": "#f24040",
        "E": "#f24040",
        "N": "#40bf40",
        "Q": "#40bf40",
        "S": "#40bf40",
        "T": "#40bf40",
        "G": "#bf40bf",
        "P": "#bf40bf",
    }


def render_logo(
    logo_df: pd.DataFrame,
    aligned_reference: str,
    mapped_ref_labels: list[str | None],
    png_path: Path,
    svg_path: Path,
    title: str,
    tick_step: int,
) -> None:
    if logomaker is None:
        raise RuntimeError(
            "logomaker is not installed. Install with: conda install -c conda-forge logomaker"
        )

    aln_len = len(aligned_reference)
    if len(mapped_ref_labels) != aln_len:
        raise ValueError("Reference-label map length does not match aligned reference length.")

    fig, (ax_logo, ax_ref) = plt.subplots(
        nrows=2,
        ncols=1,
        figsize=(max(12, aln_len * 0.28), 6.5),
        sharex=True,
        gridspec_kw={"height_ratios": [5, 1]},
    )

    logo = logomaker.Logo(logo_df, ax=ax_logo, color_scheme=_clustal_color_scheme())
    logo.style_spines(spines=["left", "bottom"], visible=True)
    logo.style_spines(spines=["top", "right"], visible=False)

    ax_logo.set_ylabel("Frequency")
    ax_logo.set_title(title)
    ax_logo.set_xlim(-0.5, aln_len - 0.5)

    ax_ref.set_ylim(0, 1)
    ax_ref.set_yticks([])
    ax_ref.set_ylabel("Ref", rotation=0, labelpad=14)

    for i, aa in enumerate(aligned_reference):
        color = "#777777" if aa == "-" else "#222222"
        ax_ref.text(i, 0.62, aa, ha="center", va="center", fontsize=7, color=color)

    step = max(1, tick_step)
    ticks = list(range(0, aln_len, step))
    tick_labels = []
    for pos in ticks:
        aln_label = str(pos + 1)
        ref_label = mapped_ref_labels[pos] if mapped_ref_labels[pos] is not None else "-"
        tick_labels.append(f"{aln_label}\n{ref_label}")

    ax_ref.set_xticks(ticks)
    ax_ref.set_xticklabels(tick_labels, fontsize=8)
    ax_ref.set_xlabel("Alignment position / Reference residue")

    for spine in ("top", "left", "right"):
        ax_ref.spines[spine].set_visible(False)

    fig.tight_layout()
    png_path.parent.mkdir(parents=True, exist_ok=True)
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(png_path), dpi=180, bbox_inches="tight")
    fig.savefig(str(svg_path), bbox_inches="tight")
    plt.close(fig)


def build_default_title(csv_path: Path, pdb_path: Path, chain_id: str) -> str:
    return f"Sequence logo: {csv_path.stem} vs {pdb_path.stem} chain {chain_id}"


def main() -> None:
    args = parse_args()

    csv_path = Path(args.csv)
    pdb_path = Path(args.pdb)
    outdir = Path(args.outdir)

    seq_df, csv_records = load_csv_sequences(
        csv_path=csv_path,
        sequence_column=args.sequence_column,
        id_column=args.id_column,
    )

    reference_sequence, reference_labels = extract_reference_sequence(
        pdb_path=pdb_path,
        chain_id=args.chain,
    )

    reference_id = f"reference|{pdb_path.stem}|chain_{args.chain}"

    if not args.allow_gaps:
        # Enforce gap-free comparison: all designed sequences must match reference length.
        assert_equal_lengths(reference_sequence, csv_records)

    combined_records: OrderedDict[str, str] = OrderedDict()
    combined_records[reference_id] = reference_sequence
    combined_records.update(csv_records)

    raw_fasta_path = outdir / f"{args.prefix}_sequences.fasta"
    aligned_fasta_path = outdir / f"{args.prefix}_aligned.fasta"
    logo_png_path = outdir / f"{args.prefix}_logo.png"
    logo_svg_path = outdir / f"{args.prefix}_logo.svg"

    write_fasta(combined_records, raw_fasta_path)

    if args.allow_gaps:
        aligned_records = star_align(combined_records)
    else:
        # In strict mode, keep sequences as-is to avoid introducing alignment gaps.
        aligned_records = combined_records
    write_fasta(aligned_records, aligned_fasta_path)

    aligned_reference = aligned_records[reference_id]
    aligned_designs = [
        aligned_seq
        for name, aligned_seq in aligned_records.items()
        if name != reference_id
    ]

    logo_df = compute_logo_matrix(aligned_designs)
    mapped_ref_labels = map_alignment_to_reference_labels(aligned_reference, reference_labels)

    title = args.title or build_default_title(csv_path, pdb_path, args.chain)
    render_logo(
        logo_df=logo_df,
        aligned_reference=aligned_reference,
        mapped_ref_labels=mapped_ref_labels,
        png_path=logo_png_path,
        svg_path=logo_svg_path,
        title=title,
        tick_step=args.tick_step,
    )

    print(f"Loaded {len(seq_df)} sequences from CSV column '{args.sequence_column}'.")
    print(f"Reference length: {len(reference_sequence)} residues (chain {args.chain}).")
    print(f"Alignment length: {len(aligned_reference)} columns.")
    print(
        "Mode: gapped alignment enabled (--allow-gaps)."
        if args.allow_gaps
        else "Mode: strict gap-free comparison (default)."
    )
    print(f"Wrote FASTA: {raw_fasta_path}")
    print(f"Wrote aligned FASTA: {aligned_fasta_path}")
    print(f"Wrote logo PNG: {logo_png_path}")
    print(f"Wrote logo SVG: {logo_svg_path}")


if __name__ == "__main__":
    main()
