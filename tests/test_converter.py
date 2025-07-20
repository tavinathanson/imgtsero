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
        self.assertIn("No serological equivalent found for molecular allele", str(cm.exception))
    
    def test_c_locus_conversion(self):
        """Test C locus conversion with Cw nomenclature."""
        # Test molecular to serological (C* to Cw)
        result = self.converter.convert("C*14:02", "s")
        self.assertEqual(result, "Cw14")
        
        result = self.converter.convert("C*12:02", "s")
        self.assertEqual(result, "Cw12")
        
        # Test with 4-field C allele
        result = self.converter.convert("C*14:02:01:01", "s")
        self.assertEqual(result, "Cw14")
        
        # Test serological to molecular (Cw to C*)
        result = self.converter.convert("Cw14", "m")
        self.assertIsInstance(result, list)
        self.assertIn("C*14:02", result)
        
        result = self.converter.convert("Cw12", "m")
        self.assertIsInstance(result, list)
        self.assertIn("C*12:02", result)
        
        # Test auto-detect for C locus
        result = self.converter.convert("C*14:02")
        self.assertEqual(result, "Cw14")
        
        result = self.converter.convert("Cw14")
        self.assertIsInstance(result, list)
        self.assertIn("C*14:02", result)
    
    def test_broad_split_conversion(self):
        """Test broad/split antigen conversion functionality."""
        # Test normal behavior (should work as before)
        result = self.converter.convert("A*02:01", "s")
        self.assertEqual(result, "A2")
        
        result = self.converter.convert("A*02:03", "s")
        self.assertEqual(result, "A203")
        
        # Test prefer_broad=True
        result = self.converter.convert("A*02:03", "s", prefer_broad=True)
        self.assertEqual(result, "A2")
        
        # Test that prefer_broad doesn't affect alleles that aren't splits
        result = self.converter.convert("A*02:01", "s", prefer_broad=True)
        self.assertEqual(result, "A2")  # A*02:01 maps to A2 (broad), so no change
        
        # Test include_broad=True for serological to molecular
        a2_normal = self.converter.convert("A2", "m")
        a2_broad = self.converter.convert("A2", "m", include_broad=True)
        
        # Should have more alleles with include_broad
        self.assertGreater(len(a2_broad), len(a2_normal))
        
        # Normal A2 alleles should be included in broad A2
        for allele in a2_normal:
            self.assertIn(allele, a2_broad)
        
        # A203 alleles should be included in broad A2
        a203_alleles = self.converter.convert("A203", "m")
        for allele in a203_alleles:
            self.assertIn(allele, a2_broad)
        
        # A210 alleles should be included in broad A2
        a210_alleles = self.converter.convert("A210", "m")
        for allele in a210_alleles:
            self.assertIn(allele, a2_broad)
    
    def test_broad_split_auto_detect(self):
        """Test broad/split functionality with auto-detection."""
        # Auto-detect with prefer_broad
        result = self.converter.convert("A*02:03", prefer_broad=True)
        self.assertEqual(result, "A2")
        
        # Auto-detect with include_broad
        a2_normal = self.converter.convert("A2")
        a2_broad = self.converter.convert("A2", include_broad=True)
        self.assertGreater(len(a2_broad), len(a2_normal))
    
    def test_broad_split_edge_cases(self):
        """Test edge cases for broad/split functionality."""
        # Test that non-split antigens aren't affected by prefer_broad
        result_normal = self.converter.convert("A*01:01", "s")
        result_broad = self.converter.convert("A*01:01", "s", prefer_broad=True)
        self.assertEqual(result_normal, result_broad)
        
        # Test that non-broad antigens aren't affected by include_broad
        a1_normal = self.converter.convert("A1", "m")
        a1_broad = self.converter.convert("A1", "m", include_broad=True)
        self.assertEqual(a1_normal, a1_broad)
    
    def test_to_2field_normalization_only(self):
        """Test that _to_2field only handles molecular normalization."""
        # Test 4-field to 2-field normalization
        result = self.converter._to_2field("A*01:01:01:01")
        self.assertEqual(result, "A*01:01")
        
        # Test 3-field to 2-field normalization  
        result = self.converter._to_2field("A*01:01:01")
        self.assertEqual(result, "A*01:01")
        
        # Test 2-field remains unchanged
        result = self.converter._to_2field("A*01:01")
        self.assertEqual(result, "A*01:01")
        
        # Test that non-molecular input is returned as-is (no conversion attempted)
        result = self.converter._to_2field("A1")
        self.assertEqual(result, "A1")  # No conversion, just returns input
        
        # Test malformed input is returned as-is
        result = self.converter._to_2field("invalid")
        self.assertEqual(result, "invalid")


if __name__ == '__main__':
    unittest.main()