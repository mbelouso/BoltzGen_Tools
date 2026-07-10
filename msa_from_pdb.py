#!/usr/bin/env python3
"""Extract protein sequences from PDB files and perform multiple sequence alignment."""

import argparse
import io
import shutil
import subprocess
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# ClustalX-style ANSI color codes by residue chemical group
_COLORS = {
    "ACFGILMVWY": "\033[33m",  # hydrophobic — yellow
    "HKR": "\033[34m",         # positive — blue
    "DE": "\033[31m",          # negative — red
    "NQST": "\033[32m",        # polar — green
    "GP": "\033[35m",          # special — magenta
}
_RESET = "\033[0m"

_CHEMICAL_GROUPS = [set(g) for g in _COLORS]

# ClustalX-style RGB colors for image output
_IMG_COLORS = {
    "A": (0.94, 0.63, 0.0), "C": (0.94, 0.63, 0.0), "F": (0.94, 0.63, 0.0),
    "I": (0.94, 0.63, 0.0), "L": (0.94, 0.63, 0.0), "M": (0.94, 0.63, 0.0),
    "V": (0.94, 0.63, 0.0), "W": (0.94, 0.63, 0.0), "Y": (0.94, 0.63, 0.0),
    "H": (0.25, 0.25, 0.95), "K": (0.25, 0.25, 0.95), "R": (0.25, 0.25, 0.95),
    "D": (0.95, 0.25, 0.25), "E": (0.95, 0.25, 0.25),
    "N": (0.25, 0.75, 0.25), "Q": (0.25, 0.75, 0.25),
    "S": (0.25, 0.75, 0.25), "T": (0.25, 0.75, 0.25),
    "G": (0.75, 0.25, 0.75), "P": (0.75, 0.25, 0.75),
    "-": (0.91, 0.91, 0.91), "X": (0.63, 0.63, 0.63),
}


def _residue_color(aa: str) -> str:
    for group, code in zip(_CHEMICAL_GROUPS, _COLORS.values()):
        if aa in group:
            return code
    return ""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract sequences from PDB files and run multiple sequence alignment."
    )
    parser.add_argument("folder", help="Path to directory containing PDB files")
    parser.add_argument("chain", help="Chain ID to extract (e.g. A)")
    parser.add_argument(
        "--output", default="msa_output.fasta",
        help="Output FASTA file for aligned sequences (default: msa_output.fasta)",
    )
    parser.add_argument("--no-color", action="store_true", help="Disable ANSI color in terminal output")
    parser.add_argument(
        "--muscle", action="store_true",
        help="Use MUSCLE for MSA instead of built-in star-alignment (requires muscle on PATH)",
    )
    parser.add_argument(
        "--reference", metavar="NAME",
        help="PDB stem to place at the top of the alignment (e.g. GLP1_Rev2)",
    )
    parser.add_argument(
        "--image", nargs="?", const="msa_output.png", default=None, metavar="FILE",
        help="Save MSA as a colored image (default filename: msa_output.png)",
    )
    return parser.parse_args()


def extract_sequence(pdb_path: str, chain_id: str) -> str:
    stem = Path(pdb_path).stem
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(stem, pdb_path)
    model = structure[0]
    if chain_id not in [c.id for c in model]:
        raise ValueError(f"Chain '{chain_id}' not found in {Path(pdb_path).name}")
    chain = model[chain_id]
    ppb = PPBuilder()
    peptides = ppb.build_peptides(chain)
    if not peptides:
        raise ValueError(f"No polypeptide residues found in chain '{chain_id}' of {Path(pdb_path).name}")
    seq = "".join(str(pp.get_sequence()) for pp in peptides)
    if "X" in seq:
        print(f"  Note: {Path(pdb_path).name} contains non-standard residues (shown as 'X')")
    return seq


def extract_all_sequences(folder: str, chain_id: str) -> dict:
    pdb_files = sorted(Path(folder).glob("*.pdb"))
    if not pdb_files:
        sys.exit(f"Error: No .pdb files found in '{folder}'")
    sequences = {}
    for pdb_path in pdb_files:
        try:
            seq = extract_sequence(str(pdb_path), chain_id)
            sequences[pdb_path.stem] = seq
        except Exception as e:
            print(f"Warning: skipping {pdb_path.name}: {e}")
    if not sequences:
        sys.exit("Error: No sequences could be extracted.")
    print(f"\nExtracted {len(sequences)} sequence(s) from '{folder}':")
    max_len = max(len(n) for n in sequences)
    for name, seq in sequences.items():
        print(f"  {name:<{max_len}}  len={len(seq)}  {seq[:30]}{'...' if len(seq) > 30 else ''}")
    return sequences


def _parse_clustal_gapped(clustal_str: str) -> tuple:
    """Extract two gapped sequence strings from a BioPython CLUSTAL-format pairwise alignment."""
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


def star_align(sequences: dict) -> dict:
    names = list(sequences.keys())
    seqs = list(sequences.values())
    n = len(seqs)

    if n == 1:
        return sequences

    aligner = _build_aligner()

    # Pick center sequence: highest total pairwise score vs all others
    scores = [sum(aligner.score(seqs[i], seqs[j]) for j in range(n) if j != i) for i in range(n)]
    center_idx = scores.index(max(scores))
    center = seqs[center_idx]
    center_len = len(center)

    # Align each sequence to center; record insertions relative to center positions
    # ins_schedule[seq_idx][center_pos] = list of inserted chars before that center position
    ins_schedule = [{} for _ in range(n)]
    match_schedule = [{} for _ in range(n)]  # match_schedule[seq_idx][center_pos] = char

    for i, seq in enumerate(seqs):
        aln = aligner.align(center, seq)[0]
        t_gapped, q_gapped = _parse_clustal_gapped(aln.format("clustal"))

        center_pos = 0
        pending_ins = []
        for t_ch, q_ch in zip(t_gapped, q_gapped):
            if t_ch == "-":
                pending_ins.append(q_ch)
            else:
                if pending_ins:
                    ins_schedule[i][center_pos] = ins_schedule[i].get(center_pos, []) + pending_ins
                    pending_ins = []
                match_schedule[i][center_pos] = q_ch
                center_pos += 1
        # trailing insertions after last center residue
        if pending_ins:
            ins_schedule[i][center_len] = ins_schedule[i].get(center_len, []) + pending_ins

    # Compute max insertion length at each center position across all sequences
    max_ins = {}
    for i in range(n):
        for pos, ins_list in ins_schedule[i].items():
            max_ins[pos] = max(max_ins.get(pos, 0), len(ins_list))

    # Reconstruct final aligned strings
    final_seqs = ["" for _ in range(n)]
    for pos in range(center_len + 1):
        # Insertion columns before this center position
        ins_width = max_ins.get(pos, 0)
        if ins_width > 0:
            for i in range(n):
                ins_chars = ins_schedule[i].get(pos, [])
                padded = "".join(ins_chars) + "-" * (ins_width - len(ins_chars))
                final_seqs[i] += padded
        # Match column (skip after last center residue)
        if pos < center_len:
            for i in range(n):
                final_seqs[i] += match_schedule[i].get(pos, "-")

    return dict(zip(names, final_seqs))


def run_muscle(sequences: dict) -> dict:
    if not shutil.which("muscle"):
        print("Warning: 'muscle' not found on PATH — falling back to built-in star-alignment.")
        return star_align(sequences)

    fasta_in = "".join(f">{name}\n{seq}\n" for name, seq in sequences.items())
    # Try MUSCLE v5 syntax first, fall back to v3/v4
    for cmd in (["muscle", "-align", "/dev/stdin", "-output", "/dev/stdout"],
                ["muscle", "-in", "/dev/stdin", "-out", "/dev/stdout"]):
        result = subprocess.run(cmd, input=fasta_in, capture_output=True, text=True)
        if result.returncode == 0:
            records = list(SeqIO.parse(io.StringIO(result.stdout), "fasta"))
            return {rec.id: str(rec.seq) for rec in records}
    print("Warning: MUSCLE failed — falling back to built-in star-alignment.")
    return star_align(sequences)


def conservation_line(aligned: dict) -> str:
    seqs = list(aligned.values())
    aln_len = len(seqs[0])
    result = []
    for i in range(aln_len):
        col = [s[i] for s in seqs if s[i] != "-"]
        if not col:
            result.append(" ")
        elif len(set(col)) == 1:
            result.append("*")
        else:
            for group in _CHEMICAL_GROUPS:
                if set(col) <= group:
                    result.append(":")
                    break
            else:
                if len(col) > len(set(col)):  # at least two identical
                    result.append(".")
                else:
                    result.append(" ")
    return "".join(result)


def format_alignment_terminal(aligned: dict, use_color: bool = True) -> str:
    names = list(aligned.keys())
    seqs = list(aligned.values())
    aln_len = len(seqs[0])
    cons = conservation_line(aligned)
    name_width = max(len(n) for n in names)
    block_width = 60
    lines = [f"\nMultiple Sequence Alignment  ({len(aligned)} sequences, {aln_len} columns)\n"]

    for start in range(0, aln_len, block_width):
        end = min(start + block_width, aln_len)
        for name, seq in zip(names, seqs):
            block = seq[start:end]
            if use_color:
                colored = "".join(
                    f"{_residue_color(c)}{c}{_RESET}" if c != "-" else c for c in block
                )
            else:
                colored = block
            lines.append(f"{name:<{name_width}}  {colored}")
        lines.append(f"{' ' * name_width}  {cons[start:end]}")
        lines.append(f"{' ' * name_width}  {start + 1:<5}   {end}")
        lines.append("")

    return "\n".join(lines)


def _find_psipred() -> str | None:
    for cmd in ("runpsipred_single", "runpsipred"):
        if shutil.which(cmd):
            return cmd
    return None


def _parse_psipred_horiz(path: str) -> str:
    parts = []
    with open(path) as f:
        for line in f:
            if line.startswith("Pred:"):
                parts.append(line[5:].strip())
    return "".join(parts)


def _parse_psipred_ss2(path: str) -> str:
    parts = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split()
            if len(cols) >= 3 and cols[2] in "HEC":
                parts.append(cols[2])
    return "".join(parts)


def _run_psipred(seq: str, name: str, cmd: str) -> str | None:
    import os
    import tempfile
    safe = "".join(c if c.isalnum() or c in "_-" else "_" for c in name)
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta = os.path.join(tmpdir, f"{safe}.fasta")
        with open(fasta, "w") as f:
            f.write(f">{safe}\n{seq}\n")
        try:
            result = subprocess.run(
                [cmd, fasta], capture_output=True, text=True,
                cwd=tmpdir, timeout=120,
            )
        except (subprocess.TimeoutExpired, OSError) as e:
            print(f"  Warning: PSIPRED error for {name}: {e}")
            return None
        if result.returncode != 0:
            return None
        for parse_fn, ext in ((_parse_psipred_horiz, ".horiz"), (_parse_psipred_ss2, ".ss2")):
            candidate = os.path.join(tmpdir, f"{safe}{ext}")
            if os.path.exists(candidate):
                ss = parse_fn(candidate)
                if ss:
                    return ss
        return None


def _map_ss_to_alignment(ss: str, aligned_seq: str) -> str:
    it = iter(ss)
    return "".join("-" if aa == "-" else next(it, "C") for aa in aligned_seq)


def extract_all_ss(sequences: dict) -> dict:
    """Run PSIPRED on each sequence string. Returns {name: H/E/C string}."""
    cmd = _find_psipred()
    if cmd is None:
        print("Warning: PSIPRED not found on PATH — secondary structure will be omitted.")
        print("  Install with: conda install -c bioconda psipred")
        return {}
    print(f"Running PSIPRED ({cmd}) on {len(sequences)} sequence(s)...")
    ss_data = {}
    for name, seq in sequences.items():
        print(f"  {name}... ", end="", flush=True)
        ss = _run_psipred(seq, name, cmd)
        if ss:
            print("done")
            ss_data[name] = ss
        else:
            print("failed")
    return ss_data


def save_fasta(aligned: dict, output_path: str) -> None:
    records = [SeqRecord(Seq(seq), id=name, description="aligned") for name, seq in aligned.items()]
    SeqIO.write(records, output_path, "fasta")
    print(f"Alignment saved to '{output_path}'")


def save_image(aligned: dict, output_path: str, ss_data: dict | None = None) -> None:
    import matplotlib.pyplot as plt

    names = list(aligned.keys())
    seqs = list(aligned.values())
    n_seqs = len(seqs)
    aln_len = len(seqs[0])
    cons = conservation_line(aligned)
    has_ss = bool(ss_data and any(ss_data.get(n) for n in names))

    cell_w, cell_h = 0.28, 0.38
    seq_spacing = 1.6 if has_ss else 1.0
    name_margin = max(len(n) for n in names) * 0.072 + 0.4

    y_cons = -1.5 if has_ss else -1.2
    ylim_bottom = -2.0 if has_ss else -1.8
    ylim_top = (n_seqs - 1) * seq_spacing + 0.9

    fig_w = aln_len * cell_w + name_margin + 2.2
    fig_h = (n_seqs * seq_spacing + 2.5) * cell_h + 0.4

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_xlim(-name_margin, aln_len)
    ax.set_ylim(ylim_bottom, ylim_top)
    ax.axis("off")

    for row, (name, seq) in enumerate(zip(names, seqs)):
        y = (n_seqs - 1 - row) * seq_spacing

        for col, aa in enumerate(seq):
            color = _IMG_COLORS.get(aa.upper(), (0.63, 0.63, 0.63))
            rect = plt.Rectangle((col, y - 0.45), 1.0, 0.9,
                                  facecolor=color, edgecolor="white", linewidth=0.5)
            ax.add_patch(rect)
            ax.text(col + 0.5, y, aa, ha="center", va="center",
                    fontsize=6.5, fontweight="bold",
                    color="white" if aa != "-" else "#AAAAAA")
        ax.text(-0.15, y, name, ha="right", va="center", fontsize=8)

        if has_ss and ss_data.get(name):
            y_ss = y - 0.70
            ax.hlines(y_ss, 0, aln_len, colors="#CCCCCC", linewidths=1.0, zorder=1)
            for col, ss_char in enumerate(ss_data[name]):
                if ss_char == "H":
                    r = plt.Rectangle((col + 0.04, y_ss - 0.22), 0.92, 0.44,
                                       facecolor="#FF8080", edgecolor="white",
                                       linewidth=0.3, zorder=2)
                    ax.add_patch(r)
                elif ss_char == "E":
                    r = plt.Rectangle((col + 0.04, y_ss - 0.15), 0.92, 0.30,
                                       facecolor="#D4A000", edgecolor="white",
                                       linewidth=0.3, zorder=2)
                    ax.add_patch(r)
            ax.text(-0.15, y_ss, "SS", ha="right", va="center",
                    fontsize=6, style="italic", color="#AAAAAA")

    # Conservation row
    ax.hlines(y_cons, 0, aln_len, colors="#EEEEEE", linewidths=0.5)
    for col, c in enumerate(cons):
        text_color = "#222222" if c == "*" else "#666666" if c == ":" else "#BBBBBB"
        ax.text(col + 0.5, y_cons, c, ha="center", va="center",
                fontsize=6.5, fontweight="bold", color=text_color)
    ax.text(-0.15, y_cons, "Conservation", ha="right", va="center",
            fontsize=7, style="italic", color="#555555")

    # Column position ticks every 5 residues
    y_ticks = (n_seqs - 1) * seq_spacing + 0.65
    for col in range(4, aln_len, 5):
        ax.text(col + 0.5, y_ticks, str(col + 1),
                ha="center", va="bottom", fontsize=5.5, color="#555555")

    # Legend
    legend_x = aln_len + 0.3
    ly = (n_seqs - 1) * seq_spacing - 0.2
    ax.text(legend_x, ly + 0.5, "Residue", ha="left", va="center",
            fontsize=7, fontweight="bold", color="#333333", clip_on=False)
    for label, color in [
        ("Hydrophobic", _IMG_COLORS["L"]), ("Positive", _IMG_COLORS["K"]),
        ("Negative",    _IMG_COLORS["D"]), ("Polar",    _IMG_COLORS["S"]),
        ("Special G/P", _IMG_COLORS["G"]), ("Gap",      _IMG_COLORS["-"]),
    ]:
        r = plt.Rectangle((legend_x, ly - 0.3), 0.7, 0.6,
                           facecolor=color, edgecolor="white", linewidth=0.5, clip_on=False)
        ax.add_patch(r)
        ax.text(legend_x + 0.85, ly, label, ha="left", va="center",
                fontsize=6.5, clip_on=False)
        ly -= 0.85

    if has_ss:
        ly -= 0.4
        ax.text(legend_x, ly, "2° Structure (PSIPRED)", ha="left", va="center",
                fontsize=7, fontweight="bold", color="#333333", clip_on=False)
        ly -= 0.75
        r = plt.Rectangle((legend_x, ly - 0.22), 0.7, 0.44,
                           facecolor="#FF8080", edgecolor="white", linewidth=0.3, clip_on=False)
        ax.add_patch(r)
        ax.text(legend_x + 0.85, ly, "Helix (H)", ha="left", va="center",
                fontsize=6.5, clip_on=False)
        ly -= 0.75
        r = plt.Rectangle((legend_x, ly - 0.15), 0.7, 0.30,
                           facecolor="#D4A000", edgecolor="white", linewidth=0.3, clip_on=False)
        ax.add_patch(r)
        ax.text(legend_x + 0.85, ly, "Strand (E)", ha="left", va="center",
                fontsize=6.5, clip_on=False)
        ly -= 0.75
        ax.hlines(ly, legend_x, legend_x + 0.7, colors="#CCCCCC", linewidths=1.0, clip_on=False)
        ax.text(legend_x + 0.85, ly, "Coil (C)", ha="left", va="center",
                fontsize=6.5, clip_on=False)

    plt.tight_layout(pad=0.2)
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Image saved to '{output_path}'")


def main() -> None:
    args = parse_args()

    folder = Path(args.folder)
    if not folder.is_dir():
        sys.exit(f"Error: '{args.folder}' is not a valid directory.")

    sequences = extract_all_sequences(str(folder), args.chain)

    # Reorder so the reference sequence appears first
    if args.reference:
        ref_key = Path(args.reference).stem  # tolerate ".pdb" suffix
        if ref_key in sequences:
            sequences = {ref_key: sequences[ref_key]} | {k: v for k, v in sequences.items() if k != ref_key}
        else:
            print(f"Warning: reference '{ref_key}' not found in extracted sequences — order unchanged.")

    if len(sequences) == 1:
        print("\nOnly 1 sequence found — no alignment needed.")
        save_fasta(sequences, args.output)
        return

    print(f"\nRunning {'MUSCLE' if args.muscle else 'star-alignment'} MSA on {len(sequences)} sequences...")
    aligned = run_muscle(sequences) if args.muscle else star_align(sequences)

    print(format_alignment_terminal(aligned, use_color=not args.no_color))
    save_fasta(aligned, args.output)

    if args.image:
        ss_raw = extract_all_ss(sequences)  # PSIPRED on ungapped sequences
        ss_data = {name: _map_ss_to_alignment(ss, aligned[name])
                   for name, ss in ss_raw.items() if ss}
        save_image(aligned, args.image, ss_data if ss_data else None)


if __name__ == "__main__":
    main()
