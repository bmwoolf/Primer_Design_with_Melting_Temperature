![Banner](assets/github_banner.png)

# Primer_Design_with_Melting_Temperature

## Challenge 2: Primer Design with Melting Temperature

### Overview

This project designs **PCR primers** from a given DNA sequence.  
The primers must:

- Be **18–24 base pairs** long
- Have **melting temperatures (Tm)** between **55–65°C**
- Optionally avoid hairpins and self-dimers

The output includes:
- Forward + reverse primer sequences
- Their Tm values
- Start/end positions on the target DNA

---

### DNA Input: pET28a::Melittin Plasmid

We use a synthetic plasmid as input:
- **Backbone**: Standard pET-28a (~5.3 kb), used for protein expression in *E. coli*
- **Insert**: The **Melittin gene**, a small toxic peptide from bee venom
- **Site**: Inserted at the **NcoI–XhoI multiple cloning site (MCS)**

This mimics a real-world cloning scenario for expressing antivenom-related genes.

---

## Files

- `pET28a_Melittin.fasta` — input plasmid sequence
- `primer_design.py` — primer selection script
- (Optional) `test_primer_design.py` — validation tests

---

## ▶Usage

```bash
python primer_design.py --fasta pET28a_Melittin.fasta --amplicon_length 500
```