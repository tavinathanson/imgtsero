"""Tests for HLA converter functionality."""

import unittest
import os
from imgtsero.converter import HLAConverter, convert, HLAConversionError
from tests.test_setup import test_data_manager


class TestHLAConverter(unittest.TestCase):
    """Test cases for HLA converter."""
    
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
        self.converter = HLAConverter(self.test_data_dir)
    
    def test_to_serological(self):
        """Test conversion to serological format with real data."""
        # Test known molecular to serological mappings
        result = self.converter.convert("A*01:01", "s")
        self.assertEqual(result, "A1")
        
        result = self.converter.convert("B*27:05", "s")
        self.assertEqual(result, "B27")
        
        # Test with 4-field allele
        result = self.converter.convert("A*01:01:01:01", "s")
        self.assertEqual(result, "A1")
    
    def test_to_molecular(self):
        """Test conversion to molecular format with real data."""
        result = self.converter.convert("A1", "m")
        self.assertIsInstance(result, list)
        self.assertGreater(len(result), 0)
        self.assertIn("A*01:01", result)
        
        # Test B27 which should have multiple molecular equivalents
        result = self.converter.convert("B27", "m")
        self.assertIsInstance(result, list)
        self.assertGreater(len(result), 0)
        self.assertIn("B*27:05", result)
        
        # Already molecular - should return 2-field format
        result = self.converter.convert("A*02:01:01:01", "m")
        self.assertEqual(result, ["A*02:01"])
    
    def test_auto_detect(self):
        """Test auto-detection and conversion to opposite format."""
        # Serological -> molecular
        result = self.converter.convert("A1")
        self.assertIsInstance(result, list)
        self.assertIn("A*01:01", result)
        
        # Molecular -> serological
        result = self.converter.convert("A*01:01")
        self.assertEqual(result, "A1")
        
        # 4-field molecular -> serological
        result = self.converter.convert("B*27:05:01:01")
        self.assertEqual(result, "B27")
    
    def test_invalid_target_format(self):
        """Test invalid target format raises error."""
        with self.assertRaises(ValueError):
            self.converter.convert("A*01:01", "invalid_format")
    
    def test_convenience_function(self):
        """Test the convenience convert function with real data."""
        # Test auto-detect
        result = convert("A1", data_dir=self.test_data_dir)
        self.assertIsInstance(result, list)
        self.assertIn("A*01:01", result)
        
        # Test explicit format
        result = convert("A*01:01", "s", self.test_data_dir)
        self.assertEqual(result, "A1")
    
    def test_real_allele_existence(self):
        """Test that conversions work with alleles that actually exist in the data."""
        # Get some real alleles from the parser
        parser = self.converter.parser
        a_alleles = list(parser.get_alleles_for_locus('A'))[:5]
        
        for allele in a_alleles:
            # Test molecular format (should normalize to 2-field)
            result_molecular = self.converter.convert(allele, "m")
            self.assertIsInstance(result_molecular, list)
            self.assertEqual(len(result_molecular), 1)
            # Should be 2-field format
            normalized = result_molecular[0]
            self.assertTrue(normalized.startswith(allele.split(':')[0] + ':' + allele.split(':')[1]))
            
            # Test that the allele format is valid
            self.assertIn('*', allele)
            self.assertIn(':', allele)
    
    def test_unrecognized_molecular_allele(self):
        """Test error handling for unrecognized molecular alleles."""
        # Test completely invalid molecular allele
        with self.assertRaises(HLAConversionError) as cm:
            self.converter.convert("A*99:99", "s")
        self.assertIn("Unrecognized molecular allele", str(cm.exception))
        
        # Test auto-detect with invalid molecular allele
        with self.assertRaises(HLAConversionError):
            self.converter.convert("B*99:99")
    
    def test_unrecognized_serological_allele(self):
        """Test error handling for unrecognized serological alleles."""
        # Test completely invalid serological allele
        with self.assertRaises(HLAConversionError) as cm:
            self.converter.convert("A999", "m")
        self.assertIn("Unrecognized serological allele", str(cm.exception))
        
        # Test auto-detect with invalid serological allele
        with self.assertRaises(HLAConversionError):
            self.converter.convert("B999")
    
    def test_invalid_format(self):
        """Test error handling for invalid HLA formats."""
        # Test completely invalid format
        with self.assertRaises(HLAConversionError) as cm:
            self.converter.convert("not-an-hla-allele")
        self.assertIn("Invalid HLA format", str(cm.exception))
        
        # Test malformed molecular format
        with self.assertRaises(HLAConversionError):
            self.converter.convert("A*invalid")
    
    def test_molecular_without_serological_equivalent(self):
        """Test molecular alleles that exist but have no serological equivalent."""
        # Use a null allele which exists in the data but has no serological equivalent
        # These alleles exist in WMDA data but map to serological "0" (no equivalent)
        with self.assertRaises(HLAConversionError) as cm:
            self.converter.convert("A*01:11N", "s")
        self.assertIn("Unrecognized molecular allele", str(cm.exception))


if __name__ == '__main__':
    unittest.main()