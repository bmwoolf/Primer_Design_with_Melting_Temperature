import argparse
import sys
from primer_design.core import load_fasta_sequence, primer_design, print_results


# Main function pipeline to design primers
def main():
    parser = argparse.ArgumentParser(description="Design PCR primers with melting temperature constraints")
    parser.add_argument("--fasta", required=True, help="Input FASTA file containing DNA sequence")
    parser.add_argument("--amplicon_length", type=int, default=500, help="Target amplicon length in base pairs (default: 500)")
    args = parser.parse_args()
    try:
        sequence = load_fasta_sequence(args.fasta)
        results = primer_design(sequence, args.amplicon_length)
        print_results(results)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main() 