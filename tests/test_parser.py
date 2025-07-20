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
    
    def test_c_locus_parsing(self):
        """Test C locus parsing with special Cw nomenclature."""
        # Test that C locus alleles are parsed
        c_alleles = self.parser.get_alleles_for_locus('C')
        self.assertGreater(len(c_alleles), 100)  # Should have many C alleles
        
        # Test C locus molecular to serological mapping
        sero = self.parser.get_molecular_to_serological('C*14:02')
        self.assertEqual(sero, 'Cw14')
        
        sero = self.parser.get_molecular_to_serological('C*12:02')
        self.assertEqual(sero, 'Cw12')
        
        # Test Cw serological to molecular mapping
        molecular_alleles = self.parser.get_serological_mapping('Cw14')
        self.assertIsInstance(molecular_alleles, list)
        self.assertIn('C*14:02', molecular_alleles)
        
        molecular_alleles = self.parser.get_serological_mapping('Cw12')
        self.assertIsInstance(molecular_alleles, list)
        self.assertIn('C*12:02', molecular_alleles)
    
    def test_broad_split_parsing(self):
        """Test broad/split antigen parsing from rel_ser_ser.txt."""
        # Test that broad/split relationships are loaded
        a2_splits = self.parser.get_split_antigens('A2')
        self.assertIsInstance(a2_splits, list)
        self.assertIn('A203', a2_splits)
        self.assertIn('A210', a2_splits)
        
        # Test reverse mapping (split to broad)
        broad = self.parser.get_broad_antigen('A203')
        self.assertEqual(broad, 'A2')
        
        broad = self.parser.get_broad_antigen('A210')
        self.assertEqual(broad, 'A2')
        
        # Test B locus broad/split relationships
        b5_splits = self.parser.get_split_antigens('B5')
        self.assertIn('B51', b5_splits)
        self.assertIn('B52', b5_splits)
        
        broad = self.parser.get_broad_antigen('B51')
        self.assertEqual(broad, 'B5')
    
    def test_broad_split_serological_mapping(self):
        """Test serological mapping with broad/split functionality."""
        # Test normal mapping
        a2_normal = self.parser.get_serological_mapping('A2', include_broad=False)
        self.assertIsInstance(a2_normal, list)
        self.assertGreater(len(a2_normal), 0)
        
        # Test with include_broad=True
        a2_broad = self.parser.get_serological_mapping('A2', include_broad=True)
        self.assertIsInstance(a2_broad, list)
        self.assertGreater(len(a2_broad), len(a2_normal))
        
        # Verify that A203 and A210 alleles are included in broad A2
        a203_alleles = self.parser.get_serological_mapping('A203')
        a210_alleles = self.parser.get_serological_mapping('A210')
        
        for allele in a203_alleles:
            self.assertIn(allele, a2_broad)
        for allele in a210_alleles:
            self.assertIn(allele, a2_broad)
    
    def test_prefer_broad_functionality(self):
        """Test prefer_broad functionality in molecular to serological mapping."""
        # Test normal behavior
        sero_normal = self.parser.get_molecular_to_serological('A*02:03', prefer_broad=False)
        self.assertEqual(sero_normal, 'A203')
        
        # Test prefer_broad=True
        sero_broad = self.parser.get_molecular_to_serological('A*02:03', prefer_broad=True)
        self.assertEqual(sero_broad, 'A2')
        
        # Test that alleles already mapping to broad aren't affected
        sero_normal = self.parser.get_molecular_to_serological('A*02:01', prefer_broad=False)
        sero_broad = self.parser.get_molecular_to_serological('A*02:01', prefer_broad=True)
        self.assertEqual(sero_normal, sero_broad)  # Both should be 'A2'


if __name__ == '__main__':
    unittest.main()