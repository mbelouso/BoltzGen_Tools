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
