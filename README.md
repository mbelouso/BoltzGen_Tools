# BoltzGen_Tools
Utility tools for BoltzGen Analysis and Generation

---

## Environment Setup

Create and activate the conda environment:

```bash
conda create -n BoltzGen_Tools python=3.11
conda activate BoltzGen_Tools
conda install -c conda-forge pandas matplotlib numpy
conda install -c conda-forge biopython
conda install -c conda-forge logomaker
```

To enable PSIPRED secondary structure prediction in the image output:

```bash
conda install -c conda-forge -c bioconda tcsh psipred
```

---

## msa_from_pdb.py

Extracts the primary sequence from every PDB file in a folder, aligns them with a multiple sequence alignment (MSA), and optionally saves a coloured image of the result.

### Basic usage

```bash
python msa_from_pdb.py <folder> <chain>
```

**Example** — align all sequences in `ranked_designs` using chain A:

```bash
python msa_from_pdb.py ./ranked_designs A
```

### Options

| Flag | Description |
|---|---|
| `--output FILE` | FASTA file for aligned sequences (default: `msa_output.fasta`) |
| `--reference NAME` | PDB stem to place at the top of the alignment (e.g. `GLP1_Rev2`) |
| `--image [FILE]` | Save a coloured PNG image of the MSA (default filename: `msa_output.png`) |
| `--muscle` | Use MUSCLE for alignment instead of the built-in star-alignment (requires `muscle` on PATH) |
| `--no-color` | Disable ANSI colour in the terminal output |

### Examples

```bash
# Save aligned FASTA and a coloured image, with a reference sequence at the top
python msa_from_pdb.py ./ranked_designs A --reference GLP1_Rev2 --image

# Custom output paths
python msa_from_pdb.py ./ranked_designs A --output my_alignment.fasta --image my_alignment.png

# Pipe plain-text alignment to a file
python msa_from_pdb.py ./ranked_designs A --no-color > alignment.txt
```

### Image output

The image shows:
- Residues coloured by chemical group (ClustalX scheme)
- A PSIPRED secondary structure strip below each sequence (helix = coral, strand = gold, coil = grey line) — requires PSIPRED to be installed
- A conservation line at the bottom
- A colour legend on the right

---

## csv_logo_from_pdb.py

Reads amino-acid sequences from a CSV column, extracts a reference sequence from a chosen PDB chain, performs a strict gap-free comparison (equal-length check), writes FASTA files, and generates a single sequence logo aligned to the reference sequence.

### Basic usage

```bash
conda run -n BoltzGen_Tools python csv_logo_from_pdb.py \
	--csv <path/to/all_designs.csv> \
	--pdb <path/to/reference.pdb> \
	--chain A
```

**Example** using this repository inputs:

```bash
conda run -n BoltzGen_Tools python csv_logo_from_pdb.py \
	--csv GLP1_binder_partialdiffusion_fampnn_4Designs_500Seq/results/all_designs.csv \
	--pdb GLP1_binder_partialdiffusion_fampnn_4Designs_500Seq/inputs/GLP1_Rev2.pdb \
	--chain A \
	--outdir GLP1_binder_partialdiffusion_fampnn_4Designs_500Seq/results \
	--prefix all_designs_refA
```

### Outputs

- `<prefix>_sequences.fasta`: unaligned FASTA containing reference + CSV sequences
- `<prefix>_aligned.fasta`: aligned FASTA (MSA)
- `<prefix>_logo.png`: sequence logo with aligned reference track
- `<prefix>_logo.svg`: vector sequence logo with aligned reference track

### Length behavior

- The script enforces that every designed sequence has the same length as the reference sequence.
- If any sequence length differs, it exits with a clear error listing mismatched entries.
- No gaps are introduced in the reference sequence when building the logo.
- Use `--allow-gaps` to enable built-in star alignment and allow sequences of different lengths.

### Key options

| Flag | Description |
|---|---|
| `--sequence-column NAME` | CSV column with sequences (default: `sequence`) |
| `--id-column NAME` | CSV column for FASTA identifiers (default: `description`) |
| `--outdir PATH` | Directory for outputs |
| `--prefix NAME` | Prefix for output filenames (default: `csv_logo`) |
| `--tick-step N` | Axis tick spacing for alignment/reference numbering (default: `5`) |
| `--title TEXT` | Optional custom logo title |
| `--allow-gaps` | Enable gapped star alignment instead of strict equal-length mode |
