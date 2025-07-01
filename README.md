![Banner](assets/github_banner.png)

# Challenge 2: Primer Design with Melting Temperature

### Overview

This project designs **PCR primers** from a given DNA sequence for protein expression in E. coli. This specific plasmid is designed for expressing the Melittin peptide in E. coli using the T7 RNA polymerase system.

The primers must:
- Be **18–24 base pairs** long
- Have **melting temperatures (Tm)** between **55–65°C**
- Optionally avoid hairpins and self-dimers

The output includes:
- Forward + reverse primer sequences
- Their Tm values
- Start/end positions on the target DNA

### DNA Input: pET28a::Melittin Plasmid

We use a synthetic plasmid as input:
- **Backbone**: Standard pET-28a (~5.3 kb), used for protein expression in *E. coli*
- **Insert**: The **Melittin gene**, a small toxic peptide from bee venom
- **Site**: Inserted at the **NcoI–XhoI multiple cloning site (MCS)**

This mimics a real-world cloning scenario for expressing antivenom-related genes.

## Files

- `pET28a_Melittin.fasta` — input plasmid sequence
- `primer_design/` — main package code
- `tests/` — unit tests

## Installation & Usage

From the project root directory:

1. **(Recommended) Create and activate a virtual environment:**
   ```bash
   python3 -m venv env
   source env/bin/activate  # On Windows use: env\Scripts\activate
   ```

2. **Install the package in editable mode:**
   ```bash
   pip install -e .
   ```

3. **Run the primer design CLI:**
   Using Python directly:
   ```bash
   python -m primer_design.cli --fasta pET28a_Melittin.fasta --amplicon_length 500
   ```
   - Replace `pET28a_Melittin.fasta` with your own FASTA file if needed.

4. **Run all unit tests:**
   ```bash
   python -m unittest discover tests
   ```

## License

MIT