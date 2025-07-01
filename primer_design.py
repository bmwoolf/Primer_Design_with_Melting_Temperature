#!/usr/bin/env python3
"""
Primer Design with Melting Temperature Calculator

This script designs PCR primers from a DNA sequence with specific constraints:
- Primer length: 18-24 bp
- Melting temperature: 55-65°C
- Target amplicon: 500 bp
"""

import argparse
import sys
from typing import Tuple, Optional
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils.MeltingTemp import Tm_NN
from Bio.SeqUtils import gc_fraction
import time
import datetime


def load_fasta_sequence(fasta_file: str) -> str:
    """
    Load DNA sequence from FASTA file.
    
    Args:
        fasta_file: Path to FASTA file
        
    Returns:
        DNA sequence as string
        
    Raises:
        FileNotFoundError: If FASTA file doesn't exist
        ValueError: If FASTA file is empty or invalid
    """
    try:
        with open(fasta_file, 'r') as handle:
            record = next(SeqIO.parse(handle, "fasta"))
            sequence = str(record.seq).upper()
            
            if not sequence:
                raise ValueError("FASTA file contains empty sequence")
                
            # Remove any whitespace or newlines
            sequence = ''.join(sequence.split())
            
            # Remove non-ATCG characters and warn if any are found
            valid_bases = set('ATCG')
            cleaned_sequence = ''.join([base for base in sequence if base in valid_bases])
            if len(cleaned_sequence) != len(sequence):
                print(f"Warning: Removed {len(sequence) - len(cleaned_sequence)} non-ATCG characters from sequence.", file=sys.stderr)
            sequence = cleaned_sequence
            
            if not sequence:
                raise ValueError("FASTA file contains no valid DNA bases after cleaning.")
            
            return sequence
            
    except FileNotFoundError:
        raise FileNotFoundError(f"FASTA file '{fasta_file}' not found")
    except StopIteration:
        raise ValueError("FASTA file is empty or invalid format")


def calculate_tm(primer: str) -> float:
    """
    Calculate melting temperature using nearest neighbor method.
    
    Args:
        primer: DNA primer sequence
        
    Returns:
        Melting temperature in Celsius
    """
    return Tm_NN(primer)


def find_valid_primers_in_window(sequence: str, window_start: int, amplicon_length: int = 500,
                                min_length: int = 18, max_length: int = 24,
                                min_tm: float = 55.0, max_tm: float = 65.0) -> tuple:
    """
    Try to find a valid forward and reverse primer pair in a given window.
    Returns (forward_primer, forward_pos, forward_tm, reverse_primer, reverse_pos, reverse_tm) or None if not found.
    """
    window_end = window_start + amplicon_length
    window_seq = sequence[window_start:window_end]
    # Forward primer: try all 18-24 bp at window start
    for fwd_len in range(min_length, max_length + 1):
        fwd_primer = window_seq[:fwd_len]
        fwd_tm = calculate_tm(fwd_primer)
        if min_tm <= fwd_tm <= max_tm:
            # Reverse primer: try all 18-24 bp at window end
            for rev_len in range(min_length, max_length + 1):
                rev_primer_template = window_seq[-rev_len:]
                rev_primer = get_reverse_complement(rev_primer_template)
                rev_tm = calculate_tm(rev_primer_template)
                if min_tm <= rev_tm <= max_tm:
                    return (fwd_primer, window_start, fwd_tm, rev_primer, window_end - rev_len, rev_tm)
    return None


def get_reverse_complement(sequence: str) -> str:
    """
    Get reverse complement of DNA sequence.
    
    Args:
        sequence: DNA sequence
        
    Returns:
        Reverse complement sequence
    """
    return str(Seq(sequence).reverse_complement())


def primer_design(fasta_file: str, amplicon_length: int = 500) -> dict:
    """
    Slide a window across the sequence and return the first valid primer pair found.
    """
    sequence = load_fasta_sequence(fasta_file)
    print(f"Loaded DNA sequence: {len(sequence)} bp")
    print(f"Sliding window search for {amplicon_length} bp amplicons...")
    start_time = time.time()
    print(f"Start time: {datetime.datetime.now().isoformat()}")
    for window_start in range(0, len(sequence) - amplicon_length + 1):
        result = find_valid_primers_in_window(sequence, window_start, amplicon_length)
        if result:
            end_time = time.time()
            print(f"End time: {datetime.datetime.now().isoformat()}")
            print(f"Elapsed time: {end_time - start_time:.2f} seconds")
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
    print(f"End time: {datetime.datetime.now().isoformat()}")
    print(f"Elapsed time: {end_time - start_time:.2f} seconds")
    raise ValueError("No suitable primer pair found in any window.")


def print_results(results: dict):
    """
    Print primer design results in a formatted way.
    
    Args:
        results: Dictionary containing primer information
    """
    print("\n" + "="*60)
    print("PRIMER DESIGN RESULTS")
    print("="*60)
    
    # Forward primer
    fwd = results['forward_primer']
    print(f"\nFORWARD PRIMER:")
    print(f"  Sequence: 5'-{fwd['sequence']}-3'")
    print(f"  Position: {fwd['position']} to {fwd['position'] + fwd['length'] - 1}")
    print(f"  Length: {fwd['length']} bp")
    print(f"  Tm: {fwd['tm']:.1f}°C")
    print(f"  GC Content: {fwd['gc_content']:.1f}%")
    
    # Reverse primer
    rev = results['reverse_primer']
    rev_template = results['reverse_primer']['template_sequence']
    print(f"\nREVERSE PRIMER:")
    print(f"  Sequence: 5'-{rev['sequence']}-3'")
    print(f"  Template: 5'-{rev_template}-3'")
    print(f"  Position: {rev['position'] - rev['length'] + 1} to {rev['position']}")
    print(f"  Length: {rev['length']} bp")
    print(f"  Tm: {rev['tm']:.1f}°C")
    print(f"  GC Content: {rev['gc_content']:.1f}%")
    
    # Amplicon information
    amp = results['amplicon']
    print(f"\nAMPLICON:")
    print(f"  Start: {amp['start']}")
    print(f"  End: {amp['end']}")
    print(f"  Length: {amp['length']} bp (target: {amp['target_length']} bp)")
    
    # Tm difference
    tm_diff = abs(fwd['tm'] - rev['tm'])
    print(f"\nPRIMER PAIR ANALYSIS:")
    print(f"  Tm difference: {tm_diff:.1f}°C")
    if tm_diff <= 5.0:
        print(f"  ✓ Primers have similar Tm (≤5°C difference)")
    else:
        print(f"  ⚠ Tm difference >5°C - may need optimization")
    
    print("\n" + "="*60)


def main():
    """Main function to run primer design from command line."""
    parser = argparse.ArgumentParser(
        description="Design PCR primers with melting temperature constraints"
    )
    parser.add_argument(
        "--fasta", 
        required=True, 
        help="Input FASTA file containing DNA sequence"
    )
    parser.add_argument(
        "--amplicon_length", 
        type=int, 
        default=500,
        help="Target amplicon length in base pairs (default: 500)"
    )
    
    args = parser.parse_args()
    
    try:
        # Design primers
        results = primer_design(args.fasta, args.amplicon_length)
        
        # Print results
        print_results(results)
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
