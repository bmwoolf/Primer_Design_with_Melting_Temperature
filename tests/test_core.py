#!/usr/bin/env python3
"""
Unit tests for primer design functionality (refactored for package structure).
"""
import unittest
import tempfile
import os
import sys
from unittest.mock import patch, MagicMock
from primer_design.core import (
    load_fasta_sequence,
    calculate_tm,
    find_valid_primers_in_window,
    get_reverse_complement,
    primer_design,
    print_results
)

class TestLoadFastaSequence(unittest.TestCase):
    def setUp(self):
        self.valid_sequence = "ATGACATCATTTAGTTGCCAGCCATGGTACAGTGAAAAGTTCTTCTCCTTTACTCATGACCACCGTAAAG"

    def test_load_valid_fasta(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(f">test_sequence\n{self.valid_sequence}\n")
            temp_file = f.name
        try:
            sequence = load_fasta_sequence(temp_file)
            self.assertEqual(sequence, self.valid_sequence)
        finally:
            os.unlink(temp_file)

    def test_load_fasta_with_whitespace(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(f">test_sequence\n{self.valid_sequence[:20]}\n{self.valid_sequence[20:40]}\n{self.valid_sequence[40:]}\n")
            temp_file = f.name
        try:
            sequence = load_fasta_sequence(temp_file)
            self.assertEqual(sequence, self.valid_sequence)
        finally:
            os.unlink(temp_file)

    def test_load_fasta_case_insensitive(self):
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
        with self.assertRaises(FileNotFoundError):
            load_fasta_sequence("nonexistent_file.fasta")

    def test_empty_fasta(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(">empty_sequence\n\n")
            temp_file = f.name
        try:
            with self.assertRaises(ValueError):
                load_fasta_sequence(temp_file)
        finally:
            os.unlink(temp_file)

    def test_invalid_dna_bases(self):
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
    def test_calculate_tm_simple_sequence(self):
        primer = "ATGCATGCATGCATGCATGC"
        tm = calculate_tm(primer)
        self.assertIsInstance(tm, float)
        self.assertGreater(tm, 0)

    def test_calculate_tm_gc_rich(self):
        primer = "GCGCGCGCGCGCGCGCGCGC"
        tm = calculate_tm(primer)
        self.assertIsInstance(tm, float)
        self.assertGreater(tm, 0)

    def test_calculate_tm_at_rich(self):
        primer = "ATATATATATATATATATAT"
        tm = calculate_tm(primer)
        self.assertIsInstance(tm, float)
        self.assertGreater(tm, 0)

class TestFindValidPrimersInWindow(unittest.TestCase):
    def setUp(self):
        self.test_sequence = "ATGACATCATTTAGTTGCCAGCCATGGTACAGTGAAAAGTTCTTCTCCTTTACTCATGACCACCGTAAAG"

    @patch('primer_design.core.calculate_tm')
    def test_find_forward_and_reverse_primer(self, mock_tm):
        mock_tm.side_effect = lambda x: 60.0 if len(x) >= 18 else 50.0
        result = find_valid_primers_in_window(self.test_sequence * 10, 0, amplicon_length=50)
        self.assertIsInstance(result, tuple)
        self.assertGreaterEqual(len(result[0]), 18)
        self.assertLessEqual(len(result[0]), 24)
        self.assertGreaterEqual(len(result[3]), 18)
        self.assertLessEqual(len(result[3]), 24)
        self.assertEqual(result[2], 60.0)
        self.assertEqual(result[5], 60.0)

    @patch('primer_design.core.calculate_tm')
    def test_no_suitable_primer(self, mock_tm):
        mock_tm.return_value = 50.0
        result = find_valid_primers_in_window(self.test_sequence * 10, 0, amplicon_length=50)
        self.assertIsNone(result)

class TestReverseComplement(unittest.TestCase):
    def test_reverse_complement_simple(self):
        sequence = "ATGC"
        rc = get_reverse_complement(sequence)
        self.assertEqual(rc, "GCAT")

    def test_reverse_complement_longer(self):
        sequence = "ATGACATCATTTAGTTGCC"
        rc = get_reverse_complement(sequence)
        expected = "GGCAACTAAATGATGTCAT"
        self.assertEqual(rc, expected)

    def test_reverse_complement_palindrome(self):
        sequence = "ATATATAT"
        rc = get_reverse_complement(sequence)
        self.assertEqual(rc, "ATATATAT")

class TestPrimerDesign(unittest.TestCase):
    def setUp(self):
        self.test_sequence = "ATGACATCATTTAGTTGCCAGCCATGGTACAGTGAAAAGTTCTTCTCCTTTACTCATGACCACCGTAAAG" * 50
        self.temp_fasta = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
        self.temp_fasta.write(f">test_sequence\n{self.test_sequence}\n")
        self.temp_fasta.close()

    def tearDown(self):
        os.unlink(self.temp_fasta.name)

    @patch('primer_design.core.calculate_tm')
    def test_primer_design_workflow(self, mock_tm):
        mock_tm.return_value = 60.0
        sequence = load_fasta_sequence(self.temp_fasta.name)
        results = primer_design(sequence, amplicon_length=500)
        self.assertIn('forward_primer', results)
        self.assertIn('reverse_primer', results)
        self.assertIn('amplicon', results)
        fwd = results['forward_primer']
        self.assertIn('sequence', fwd)
        self.assertIn('position', fwd)
        self.assertIn('length', fwd)
        self.assertIn('tm', fwd)
        self.assertIn('gc_content', fwd)
        rev = results['reverse_primer']
        self.assertIn('sequence', rev)
        self.assertIn('template_sequence', rev)
        self.assertIn('position', rev)
        self.assertIn('length', rev)
        self.assertIn('tm', rev)
        self.assertIn('gc_content', rev)
        amp = results['amplicon']
        self.assertIn('start', amp)
        self.assertIn('end', amp)
        self.assertIn('length', amp)
        self.assertIn('target_length', amp)
        self.assertGreaterEqual(fwd['length'], 18)
        self.assertLessEqual(fwd['length'], 24)
        self.assertGreaterEqual(rev['length'], 18)
        self.assertLessEqual(rev['length'], 24)
        self.assertEqual(fwd['tm'], 60.0)
        self.assertEqual(rev['tm'], 60.0)

class TestPrintResults(unittest.TestCase):
    def setUp(self):
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
        print_results(self.test_results)
        self.assertGreater(mock_print.call_count, 10)
        print_calls = [call[0][0] for call in mock_print.call_args_list]
        printed_text = '\n'.join(print_calls)
        self.assertIn('PRIMER DESIGN RESULTS', printed_text)
        self.assertIn('FORWARD PRIMER:', printed_text)
        self.assertIn('REVERSE PRIMER:', printed_text)
        self.assertIn('AMPLICON:', printed_text)
        self.assertIn('ATGACATCATTTAGTTGCC', printed_text)
        self.assertIn('GCATGCATGCATGCATGCAT', printed_text)

class TestIntegration(unittest.TestCase):
    def setUp(self):
        self.test_sequence = "ATGACATCATTTAGTTGCCAGCCATGGTACAGTGAAAAGTTCTTCTCCTTTACTCATGACCACCGTAAAG" * 100
        self.temp_fasta = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
        self.temp_fasta.write(f">test_sequence\n{self.test_sequence}\n")
        self.temp_fasta.close()

    def tearDown(self):
        os.unlink(self.temp_fasta.name)

    def test_integration_with_real_tm_calculation(self):
        try:
            sequence = load_fasta_sequence(self.temp_fasta.name)
            results = primer_design(sequence, amplicon_length=500)
            self.assertIsInstance(results, dict)
            self.assertIn('forward_primer', results)
            self.assertIn('reverse_primer', results)
            self.assertIn('amplicon', results)
            fwd = results['forward_primer']
            rev = results['reverse_primer']
            self.assertGreaterEqual(fwd['length'], 18)
            self.assertLessEqual(fwd['length'], 24)
            self.assertGreaterEqual(rev['length'], 18)
            self.assertLessEqual(rev['length'], 24)
            self.assertGreater(fwd['tm'], 50)
            self.assertLess(fwd['tm'], 70)
            self.assertGreater(rev['tm'], 50)
            self.assertLess(rev['tm'], 70)
        except Exception as e:
            self.fail(f"Integration test failed with exception: {e}")

if __name__ == '__main__':
    unittest.main(verbosity=2) 