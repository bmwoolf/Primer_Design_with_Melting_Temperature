import sys
from typing import Tuple
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils.MeltingTemp import Tm_NN
from Bio.SeqUtils import gc_fraction
import time


# Load the FASTA sequence from the file
def load_fasta_sequence(fasta_file: str) -> str:
    try:
        with open(fasta_file, 'r') as handle:
            record = next(SeqIO.parse(handle, "fasta"))
            sequence = str(record.seq).upper()
            if not sequence:
                raise ValueError("FASTA file contains empty sequence")
            sequence = ''.join(sequence.split())
            valid_bases = set('ATCG')
            invalid_bases = set(sequence) - valid_bases
            if invalid_bases:
                raise ValueError(f"FASTA file contains invalid DNA bases: {''.join(sorted(invalid_bases))}")
            return sequence
    except FileNotFoundError:
        raise FileNotFoundError(f"FASTA file '{fasta_file}' not found")
    except StopIteration:
        raise ValueError("FASTA file is empty or invalid format")


# Calculate the melting temperature of a primer
def calculate_tm(primer: str) -> float:
    return Tm_NN(primer)


# Get the reverse complement of a sequence
def get_reverse_complement(sequence: str) -> str:
    return str(Seq(sequence).reverse_complement())


# Find valid primers in a window
def find_valid_primers_in_window(sequence: str, window_start: int, amplicon_length: int = 500,
                                min_length: int = 18, max_length: int = 24,
                                min_tm: float = 55.0, max_tm: float = 65.0) -> tuple:
    window_end = window_start + amplicon_length
    window_seq = sequence[window_start:window_end]
    for fwd_len in range(min_length, max_length + 1):
        fwd_primer = window_seq[:fwd_len]
        fwd_tm = calculate_tm(fwd_primer)
        if min_tm <= fwd_tm <= max_tm:
            for rev_len in range(min_length, max_length + 1):
                rev_primer_template = window_seq[-rev_len:]
                rev_primer = get_reverse_complement(rev_primer_template)
                rev_tm = calculate_tm(rev_primer_template)
                if min_tm <= rev_tm <= max_tm:
                    return (fwd_primer, window_start, fwd_tm, rev_primer, window_end - rev_len, rev_tm)
    return None


# Main function to design primers
def primer_design(sequence: str, amplicon_length: int = 500) -> dict:
    start_time = time.time()
    for window_start in range(0, len(sequence) - amplicon_length + 1):
        result = find_valid_primers_in_window(sequence, window_start, amplicon_length)
        if result:
            end_time = time.time()
            (fwd_primer, fwd_pos, fwd_tm, rev_primer, rev_pos, rev_tm) = result
            return {
                'forward_primer': {
                    'sequence': fwd_primer,
                    'position': fwd_pos,
                    'length': len(fwd_primer),
                    'tm': fwd_tm,
                    'gc_content': gc_fraction(fwd_primer) * 100
                },
                'reverse_primer': {
                    'sequence': rev_primer,
                    'template_sequence': sequence[rev_pos:rev_pos+len(rev_primer)],
                    'position': rev_pos,
                    'length': len(rev_primer),
                    'tm': rev_tm,
                    'gc_content': gc_fraction(sequence[rev_pos:rev_pos+len(rev_primer)]) * 100
                },
                'amplicon': {
                    'start': fwd_pos,
                    'end': rev_pos + len(rev_primer),
                    'length': (rev_pos + len(rev_primer)) - fwd_pos,
                    'target_length': amplicon_length
                },
                'search': {
                    'start_time': start_time,
                    'end_time': end_time,
                    'elapsed': end_time - start_time,
                    'window_start': window_start
                }
            }
    end_time = time.time()
    raise ValueError("No suitable primer pair found in any window.")


# Print the results in a readable format
def print_results(results: dict):
    print("\n" + "="*60)
    print("PRIMER DESIGN RESULTS")
    print("="*60)
    fwd = results['forward_primer']
    print(f"\nFORWARD PRIMER:")
    print(f"  Sequence: 5'-{fwd['sequence']}-3'")
    print(f"  Position: {fwd['position']} to {fwd['position'] + fwd['length'] - 1}")
    print(f"  Length: {fwd['length']} bp")
    print(f"  Tm: {fwd['tm']:.1f}°C")
    print(f"  GC Content: {fwd['gc_content']:.1f}%")
    rev = results['reverse_primer']
    rev_template = results['reverse_primer']['template_sequence']
    print(f"\nREVERSE PRIMER:")
    print(f"  Sequence: 5'-{rev['sequence']}-3'")
    print(f"  Template: 5'-{rev_template}-3'")
    print(f"  Position: {rev['position']} to {rev['position'] + rev['length'] - 1}")
    print(f"  Length: {rev['length']} bp")
    print(f"  Tm: {rev['tm']:.1f}°C")
    print(f"  GC Content: {rev['gc_content']:.1f}%")
    amp = results['amplicon']
    print(f"\nAMPLICON:")
    print(f"  Start: {amp['start']}")
    print(f"  End: {amp['end']}")
    print(f"  Length: {amp['length']} bp (target: {amp['target_length']} bp)")
    tm_diff = abs(fwd['tm'] - rev['tm'])
    print(f"\nPRIMER PAIR ANALYSIS:")
    print(f"  Tm difference: {tm_diff:.1f}°C")
    if tm_diff <= 5.0:
        print(f"  ✓ Primers have similar Tm (≤5°C difference)")
    else:
        print(f"  ⚠ Tm difference >5°C - may need optimization")
    print("\n" + "="*60) 