"""Tests for HLA parser functionality."""

import unittest
import os
from imgtsero.parser import HLAParser
from tests.test_setup import test_data_manager


class TestHLAParser(unittest.TestCase):
    """Test cases for HLA parser."""
    
    @classmethod
    def setUpClass(cls):
        """Set up test data once for all tests."""
        cls.test_data_dir = test_data_manager.setup_test_data()
    
    @classmethod
    def tearDownClass(cls):
        """Clean up test data."""
        test_data_manager.cleanup()
    
    def setUp(self):
        """Set up test fixtures."""
        self.parser = HLAParser(self.test_data_dir)
    
    def test_get_loci(self):
        """Test getting all available loci with real data."""
        loci = self.parser.get_loci()
        # Real IMGT data should have these standard loci
        expected_loci = ['A', 'B', 'C', 'DRB1', 'DQB1', 'DPB1']
        for locus in expected_loci[:4]:  # Test at least the main classical loci
            self.assertIn(locus, loci)
        self.assertGreater(len(loci), 5)  # Should have many loci
    
    def test_get_alleles_for_locus(self):
        """Test getting alleles for specific locus with real data."""
        a_alleles = self.parser.get_alleles_for_locus('A')
        self.assertGreater(len(a_alleles), 100)  # Real data should have many A alleles
        
        # Check that alleles have proper format
        for allele in list(a_alleles)[:5]:  # Check first 5
            self.assertTrue(allele.startswith('A*'))
            self.assertIn(':', allele)
        
        # Test non-existent locus
        empty_alleles = self.parser.get_alleles_for_locus('NONEXISTENT')
        self.assertEqual(len(empty_alleles), 0)
    
    def test_find_alleles_by_pattern(self):
        """Test pattern-based allele search with real data."""
        # Find all A*01 alleles
        results = self.parser.find_alleles_by_pattern(r'A\*01')
        self.assertGreater(len(results), 5)  # Should be multiple A*01 variants
        if results:  # Only test if we found results
            for allele in results[:5]:  # Check first few
                self.assertTrue(allele.startswith('A*01'), f"Allele {allele} doesn't start with A*01")
        
        # Find A*02:01 alleles (common allele)
        results = self.parser.find_alleles_by_pattern(r'A\*02:01')
        self.assertGreater(len(results), 0)
        # Should find A*02:01 variants
        found_common = any('A*02:01' in allele for allele in results)
        self.assertTrue(found_common)
    
    def test_serological_mapping(self):
        """Test serological to molecular mapping."""
        # Test known mapping
        molecular_alleles = self.parser.get_serological_mapping('A1')
        self.assertIsInstance(molecular_alleles, list)
        self.assertGreater(len(molecular_alleles), 0)
        
        # Test reverse mapping
        sero = self.parser.get_molecular_to_serological('A*01:01')
        self.assertEqual(sero, 'A1')
        
        # Test unknown serological type
        unknown = self.parser.get_serological_mapping('UNKNOWN99')
        self.assertEqual(unknown, [])


if __name__ == '__main__':
    unittest.main()