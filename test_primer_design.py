#!/usr/bin/env python3
"""
Unit tests for primer design functionality.

Tests cover:
- FASTA file loading
- Melting temperature calculation
- Primer finding algorithms
- Reverse complement generation
- Complete primer design workflow
"""

import unittest
import tempfile
import os
import sys
from unittest.mock import patch, MagicMock

# Import the functions to test
from primer_design import (
    load_fasta_sequence,
    calculate_tm,
    find_optimal_primer,
    get_reverse_complement,
    primer_design,
    print_results
)


class TestLoadFastaSequence(unittest.TestCase):
    """Test FASTA file loading functionality."""
    
    def setUp(self):
        """Set up test data."""
        self.valid_sequence = "ATGACATCATTTAGTTGCCAGCCATGGTACAGTGAAAAGTTCTTCTCCTTTACTCATGACCACCGTAAAG"
        
    def test_load_valid_fasta(self):
        """Test loading a valid FASTA file."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(f">test_sequence\n{self.valid_sequence}\n")
            temp_file = f.name
        
        try:
            sequence = load_fasta_sequence(temp_file)
            self.assertEqual(sequence, self.valid_sequence)
        finally:
            os.unlink(temp_file)
    
    def test_load_fasta_with_whitespace(self):
        """Test loading FASTA file with whitespace."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(f">test_sequence\n{self.valid_sequence[:20]}\n{self.valid_sequence[20:40]}\n{self.valid_sequence[40:]}\n")
            temp_file = f.name
        
        try:
            sequence = load_fasta_sequence(temp_file)
            self.assertEqual(sequence, self.valid_sequence)
        finally:
            os.unlink(temp_file)
    
    def test_load_fasta_case_insensitive(self):
        """Test loading FASTA file with mixed case."""
        mixed_case = "atgacatcattTAGTTGCCAGCCATGGTACAGTGAAAAGTTCTTCTCCTTTACTCATGACCACCGTAAAG"
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(f">test_sequence\n{mixed_case}\n")
            temp_file = f.name
        
        try:
            sequence = load_fasta_sequence(temp_file)
            self.assertEqual(sequence, mixed_case.upper())
        finally:
            os.unlink(temp_file)
    
    def test_file_not_found(self):
        """Test handling of non-existent file."""
        with self.assertRaises(FileNotFoundError):
            load_fasta_sequence("nonexistent_file.fasta")
    
    def test_empty_fasta(self):
        """Test handling of empty FASTA file."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(">empty_sequence\n\n")
            temp_file = f.name
        
        try:
            with self.assertRaises(ValueError):
                load_fasta_sequence(temp_file)
        finally:
            os.unlink(temp_file)
    
    def test_invalid_dna_bases(self):
        """Test handling of invalid DNA bases."""
        invalid_sequence = "ATGACATCATTTAGTTGCCAGCCATGGTACAGTGAAAAGTTCTTCTCCTTTACTCATGACCACCGTAAAGX"
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(f">test_sequence\n{invalid_sequence}\n")
            temp_file = f.name
        
        try:
            with self.assertRaises(ValueError):
                load_fasta_sequence(temp_file)
        finally:
            os.unlink(temp_file)


class TestCalculateTm(unittest.TestCase):
    """Test melting temperature calculation."""
    
    def test_calculate_tm_simple_sequence(self):
        """Test Tm calculation for a simple sequence."""
        primer = "ATGCATGCATGCATGCATGC"
        tm = calculate_tm(primer)
        self.assertIsInstance(tm, float)
        self.assertGreater(tm, 0)
    
    def test_calculate_tm_gc_rich(self):
        """Test Tm calculation for GC-rich sequence."""
        primer = "GCGCGCGCGCGCGCGCGCGC"
        tm = calculate_tm(primer)
        self.assertIsInstance(tm, float)
        self.assertGreater(tm, 0)
    
    def test_calculate_tm_at_rich(self):
        """Test Tm calculation for AT-rich sequence."""
        primer = "ATATATATATATATATATAT"
        tm = calculate_tm(primer)
        self.assertIsInstance(tm, float)
        self.assertGreater(tm, 0)


class TestFindOptimalPrimer(unittest.TestCase):
    """Test optimal primer finding algorithm."""
    
    def setUp(self):
        """Set up test sequence."""
        self.test_sequence = "ATGACATCATTTAGTTGCCAGCCATGGTACAGTGAAAAGTTCTTCTCCTTTACTCATGACCACCGTAAAG"
    
    @patch('primer_design.calculate_tm')
    def test_find_forward_primer(self, mock_tm):
        """Test finding forward primer."""
        # Mock Tm calculation to return values in range
        mock_tm.side_effect = lambda x: 60.0 if len(x) >= 18 else 50.0
        
        primer, pos, tm = find_optimal_primer(self.test_sequence, 0, is_forward=True)
        
        self.assertIsInstance(primer, str)
        self.assertGreaterEqual(len(primer), 18)
        self.assertLessEqual(len(primer), 24)
        self.assertEqual(tm, 60.0)
    
    @patch('primer_design.calculate_tm')
    def test_find_reverse_primer(self, mock_tm):
        """Test finding reverse primer."""
        # Mock Tm calculation to return values in range
        mock_tm.side_effect = lambda x: 60.0 if len(x) >= 18 else 50.0
        
        primer, pos, tm = find_optimal_primer(self.test_sequence, len(self.test_sequence), is_forward=False)
        
        self.assertIsInstance(primer, str)
        self.assertGreaterEqual(len(primer), 18)
        self.assertLessEqual(len(primer), 24)
        self.assertEqual(tm, 60.0)
    
    @patch('primer_design.calculate_tm')
    def test_no_suitable_primer(self, mock_tm):
        """Test handling when no suitable primer is found."""
        # Mock Tm calculation to return values outside range
        mock_tm.return_value = 50.0
        
        with self.assertRaises(ValueError):
            find_optimal_primer(self.test_sequence, 0, is_forward=True)


class TestReverseComplement(unittest.TestCase):
    """Test reverse complement generation."""
    
    def test_reverse_complement_simple(self):
        """Test reverse complement of simple sequence."""
        sequence = "ATGC"
        rc = get_reverse_complement(sequence)
        self.assertEqual(rc, "GCAT")
    
    def test_reverse_complement_longer(self):
        """Test reverse complement of longer sequence."""
        sequence = "ATGACATCATTTAGTTGCC"
        rc = get_reverse_complement(sequence)
        expected = "GGCAACTAAATGATGTCAT"
        self.assertEqual(rc, expected)
    
    def test_reverse_complement_palindrome(self):
        """Test reverse complement of palindrome."""
        sequence = "ATATATAT"
        rc = get_reverse_complement(sequence)
        self.assertEqual(rc, "ATATATAT")


class TestPrimerDesign(unittest.TestCase):
    """Test complete primer design workflow."""
    
    def setUp(self):
        """Set up test FASTA file."""
        self.test_sequence = "ATGACATCATTTAGTTGCCAGCCATGGTACAGTGAAAAGTTCTTCTCCTTTACTCATGACCACCGTAAAG" * 50
        self.temp_fasta = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
        self.temp_fasta.write(f">test_sequence\n{self.test_sequence}\n")
        self.temp_fasta.close()
    
    def tearDown(self):
        """Clean up test file."""
        os.unlink(self.temp_fasta.name)
    
    @patch('primer_design.calculate_tm')
    def test_primer_design_workflow(self, mock_tm):
        """Test complete primer design workflow."""
        # Mock Tm calculation to return values in range
        mock_tm.return_value = 60.0
        
        results = primer_design(self.temp_fasta.name, amplicon_length=500)
        
        # Check structure of results
        self.assertIn('forward_primer', results)
        self.assertIn('reverse_primer', results)
        self.assertIn('amplicon', results)
        
        # Check forward primer
        fwd = results['forward_primer']
        self.assertIn('sequence', fwd)
        self.assertIn('position', fwd)
        self.assertIn('length', fwd)
        self.assertIn('tm', fwd)
        self.assertIn('gc_content', fwd)
        
        # Check reverse primer
        rev = results['reverse_primer']
        self.assertIn('sequence', rev)
        self.assertIn('template_sequence', rev)
        self.assertIn('position', rev)
        self.assertIn('length', rev)
        self.assertIn('tm', rev)
        self.assertIn('gc_content', rev)
        
        # Check amplicon
        amp = results['amplicon']
        self.assertIn('start', amp)
        self.assertIn('end', amp)
        self.assertIn('length', amp)
        self.assertIn('target_length', amp)
        
        # Check primer constraints
        self.assertGreaterEqual(fwd['length'], 18)
        self.assertLessEqual(fwd['length'], 24)
        self.assertGreaterEqual(rev['length'], 18)
        self.assertLessEqual(rev['length'], 24)
        self.assertEqual(fwd['tm'], 60.0)
        self.assertEqual(rev['tm'], 60.0)


class TestPrintResults(unittest.TestCase):
    """Test results printing functionality."""
    
    def setUp(self):
        """Set up test results data."""
        self.test_results = {
            'forward_primer': {
                'sequence': 'ATGACATCATTTAGTTGCC',
                'position': 1000,
                'length': 19,
                'tm': 60.5,
                'gc_content': 42.1
            },
            'reverse_primer': {
                'sequence': 'GCATGCATGCATGCATGCAT',
                'template_sequence': 'ATGCATGCATGCATGCATGC',
                'position': 1500,
                'length': 20,
                'tm': 61.2,
                'gc_content': 50.0
            },
            'amplicon': {
                'start': 1000,
                'end': 1520,
                'length': 520,
                'target_length': 500
            }
        }
    
    @patch('builtins.print')
    def test_print_results(self, mock_print):
        """Test that print_results calls print with expected content."""
        print_results(self.test_results)
        
        # Check that print was called multiple times
        self.assertGreater(mock_print.call_count, 10)
        
        # Check for specific content in print calls
        print_calls = [call[0][0] for call in mock_print.call_args_list]
        printed_text = '\n'.join(print_calls)
        
        # Check for key information
        self.assertIn('PRIMER DESIGN RESULTS', printed_text)
        self.assertIn('FORWARD PRIMER:', printed_text)
        self.assertIn('REVERSE PRIMER:', printed_text)
        self.assertIn('AMPLICON:', printed_text)
        self.assertIn('ATGACATCATTTAGTTGCC', printed_text)
        self.assertIn('GCATGCATGCATGCATGCAT', printed_text)


class TestIntegration(unittest.TestCase):
    """Integration tests for the complete system."""
    
    def setUp(self):
        """Set up test FASTA file with realistic sequence."""
        # Create a longer sequence for realistic testing
        self.test_sequence = "ATGACATCATTTAGTTGCCAGCCATGGTACAGTGAAAAGTTCTTCTCCTTTACTCATGACCACCGTAAAG" * 100
        self.temp_fasta = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
        self.temp_fasta.write(f">test_sequence\n{self.test_sequence}\n")
        self.temp_fasta.close()
    
    def tearDown(self):
        """Clean up test file."""
        os.unlink(self.temp_fasta.name)
    
    def test_integration_with_real_tm_calculation(self):
        """Test integration with real Tm calculation."""
        try:
            results = primer_design(self.temp_fasta.name, amplicon_length=500)
            
            # Verify results structure
            self.assertIsInstance(results, dict)
            self.assertIn('forward_primer', results)
            self.assertIn('reverse_primer', results)
            self.assertIn('amplicon', results)
            
            # Verify primer properties
            fwd = results['forward_primer']
            rev = results['reverse_primer']
            
            self.assertGreaterEqual(fwd['length'], 18)
            self.assertLessEqual(fwd['length'], 24)
            self.assertGreaterEqual(rev['length'], 18)
            self.assertLessEqual(rev['length'], 24)
            
            # Verify Tm values are reasonable
            self.assertGreater(fwd['tm'], 50)
            self.assertLess(fwd['tm'], 70)
            self.assertGreater(rev['tm'], 50)
            self.assertLess(rev['tm'], 70)
            
        except Exception as e:
            self.fail(f"Integration test failed with exception: {e}")


if __name__ == '__main__':
    # Run tests
    unittest.main(verbosity=2)
