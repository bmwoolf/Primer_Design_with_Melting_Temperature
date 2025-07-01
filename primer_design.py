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


def find_optimal_primer(sequence: str, start_pos: int, is_forward: bool = True, 
                       min_length: int = 18, max_length: int = 24,
                       min_tm: float = 55.0, max_tm: float = 65.0) -> Tuple[str, int, float]:
    """
    Find optimal primer within specified constraints.
    
    Args:
        sequence: DNA sequence to search
        start_pos: Starting position for primer search
        is_forward: True for forward primer, False for reverse primer
        min_length: Minimum primer length
        max_length: Maximum primer length
        min_tm: Minimum melting temperature
        max_tm: Maximum melting temperature
        
    Returns:
        Tuple of (primer_sequence, position, melting_temperature)
        
    Raises:
        ValueError: If no suitable primer found
    """
    best_primer = None
    best_tm = 0.0
    best_pos = 0
    
    # For forward primer, search forward from start_pos
    # For reverse primer, search backward from start_pos
    if is_forward:
        search_range = range(start_pos, min(start_pos + 100, len(sequence) - min_length))
    else:
        search_range = range(max(0, start_pos - 100), start_pos - min_length, -1)
    
    for pos in search_range:
        for length in range(min_length, max_length + 1):
            if is_forward:
                if pos + length > len(sequence):
                    continue
                primer = sequence[pos:pos + length]
            else:
                if pos - length < 0:
                    continue
                primer = sequence[pos - length:pos]
            
            # Calculate melting temperature
            try:
                tm = calculate_tm(primer)
                
                # Check if temperature is within range
                if min_tm <= tm <= max_tm:
                    # Prefer primers closer to target temperature (60°C)
                    target_tm = 60.0
                    if best_primer is None or abs(tm - target_tm) < abs(best_tm - target_tm):
                        best_primer = primer
                        best_tm = tm
                        best_pos = pos
                        
            except Exception:
                # Skip primers that can't be calculated
                continue
    
    if best_primer is None:
        raise ValueError(f"No suitable {'forward' if is_forward else 'reverse'} primer found")
    
    return best_primer, best_pos, best_tm


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
    Design forward and reverse primers for PCR amplification.
    
    Args:
        fasta_file: Path to FASTA file containing DNA sequence
        amplicon_length: Target amplicon length in base pairs
        
    Returns:
        Dictionary containing primer information
    """
    # Load DNA sequence
    sequence = load_fasta_sequence(fasta_file)
    print(f"Loaded DNA sequence: {len(sequence)} bp")
    
    # Define amplicon region (start and end positions)
    # For this example, we'll start from position 1000 to avoid plasmid backbone
    amplicon_start = 1000
    amplicon_end = amplicon_start + amplicon_length
    
    if amplicon_end > len(sequence):
        raise ValueError(f"Amplicon end position ({amplicon_end}) exceeds sequence length ({len(sequence)})")
    
    print(f"Target amplicon: positions {amplicon_start} to {amplicon_end} ({amplicon_length} bp)")
    
    # Design forward primer (from start of amplicon)
    print("\nDesigning forward primer...")
    forward_primer, forward_pos, forward_tm = find_optimal_primer(
        sequence, amplicon_start, is_forward=True
    )
    
    # Design reverse primer (from end of amplicon, using reverse complement)
    print("Designing reverse primer...")
    reverse_primer, reverse_pos, reverse_tm = find_optimal_primer(
        sequence, amplicon_end, is_forward=False
    )
    
    # Get reverse complement of reverse primer for actual PCR use
    reverse_primer_rc = get_reverse_complement(reverse_primer)
    
    # Calculate actual amplicon length
    actual_amplicon_length = reverse_pos - forward_pos + len(forward_primer)
    
    # Prepare results
    results = {
        'forward_primer': {
            'sequence': forward_primer,
            'position': forward_pos,
            'length': len(forward_primer),
            'tm': forward_tm,
            'gc_content': gc_fraction(forward_primer) * 100
        },
        'reverse_primer': {
            'sequence': reverse_primer_rc,  # Reverse complement for PCR
            'template_sequence': reverse_primer,  # Original sequence on template
            'position': reverse_pos,
            'length': len(reverse_primer),
            'tm': reverse_tm,
            'gc_content': gc_fraction(reverse_primer) * 100
        },
        'amplicon': {
            'start': forward_pos,
            'end': reverse_pos + len(reverse_primer),
            'length': actual_amplicon_length,
            'target_length': amplicon_length
        }
    }
    
    return results


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
