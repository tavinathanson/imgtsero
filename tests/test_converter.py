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
        self.converter = HLAConverter(3610, self.test_data_dir)
    
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
        result = convert("A1", 3610, data_dir=self.test_data_dir)
        self.assertIsInstance(result, list)
        self.assertIn("A*01:01", result)
        
        # Test explicit format
        result = convert("A*01:01", 3610, "s", self.test_data_dir)
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
    
    def test_c_locus_comprehensive_bidirectional(self):
        """Test comprehensive bidirectional conversion for all C locus antigens."""
        # Test cases based on WMDA data and HLA nomenclature standards
        # According to research, ALL C locus serological antigens use Cw prefix
        test_cases = [
            ("C*01:02", "Cw1"),  # Serological 1 -> Cw1
            ("C*02:02", "Cw2"),  # Serological 2 -> Cw2  
            ("C*03:02", "Cw10"), # Serological 10 -> Cw10
            ("C*03:03", "Cw9"),  # Serological 9 -> Cw9
            ("C*03:07", "Cw3"),  # Serological 3 -> Cw3
            ("C*04:01", "Cw4"),  # Serological 4 -> Cw4
            ("C*05:01", "Cw5"),  # Serological 5 -> Cw5
            ("C*06:02", "Cw6"),  # Serological 6 -> Cw6
            ("C*07:01", "Cw7"),  # Serological 7 -> Cw7
            ("C*08:01", "Cw8"),  # Serological 8 -> Cw8
            ("C*12:02", "Cw12"), # Question mark cases -> use allele number
            ("C*14:02", "Cw14"), # Question mark cases -> use allele number
            ("C*15:02", "Cw15"), # Question mark cases -> use allele number
            ("C*16:01", "Cw16"), # Question mark cases -> use allele number
            ("C*17:01", "Cw17"), # Question mark cases -> use allele number
            ("C*18:01", "Cw18"), # Question mark cases -> use allele number
        ]
        
        for molecular, expected_sero in test_cases:
            with self.subTest(molecular=molecular, expected_sero=expected_sero):
                # Test molecular to serological
                sero_result = self.converter.convert(molecular, "s")
                self.assertEqual(sero_result, expected_sero, 
                    f"Expected {molecular} -> {expected_sero}, got {sero_result}")
                
                # Test auto-detect molecular to serological  
                auto_sero_result = self.converter.convert(molecular)
                self.assertEqual(auto_sero_result, expected_sero,
                    f"Auto-detect: Expected {molecular} -> {expected_sero}, got {auto_sero_result}")
                
                # Test serological to molecular (bidirectional)
                molecular_result = self.converter.convert(expected_sero, "m")
                self.assertIsInstance(molecular_result, list,
                    f"Expected list for {expected_sero} -> molecular, got {type(molecular_result)}")
                self.assertIn(molecular, molecular_result,
                    f"Expected {molecular} in results for {expected_sero}, got {molecular_result}")
                
                # Test auto-detect serological to molecular
                auto_molecular_result = self.converter.convert(expected_sero)
                self.assertIsInstance(auto_molecular_result, list,
                    f"Auto-detect: Expected list for {expected_sero} -> molecular, got {type(auto_molecular_result)}")
                self.assertIn(molecular, auto_molecular_result,
                    f"Auto-detect: Expected {molecular} in results for {expected_sero}, got {auto_molecular_result}")
    
    def test_c_locus_invalid_old_nomenclature(self):
        """Test that old C nomenclature without 'w' is not accepted."""
        # Test that C1, C2, etc. (without 'w') are rejected
        invalid_old_names = ["C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"]
        
        for invalid_name in invalid_old_names:
            with self.subTest(invalid_name=invalid_name):
                with self.assertRaises(HLAConversionError) as cm:
                    self.converter.convert(invalid_name, "m")
                self.assertIn("Unrecognized serological allele", str(cm.exception))
                
                # Test auto-detect also rejects it
                with self.assertRaises(HLAConversionError):
                    self.converter.convert(invalid_name)
    
    def test_broad_split_conversion(self):
        """Test broad/split antigen conversion functionality."""
        # Test normal behavior (should work as before)
        result = self.converter.convert("A*02:01", "s")
        self.assertEqual(result, "A2")
        
        result = self.converter.convert("A*02:03", "s")
        self.assertEqual(result, "A203")
        
        # Test return_broad=True
        result = self.converter.convert("A*02:03", "s", return_broad=True)
        self.assertEqual(result, "A2")
        
        # Test that return_broad doesn't affect alleles that aren't splits
        result = self.converter.convert("A*02:01", "s", return_broad=True)
        self.assertEqual(result, "A2")  # A*02:01 maps to A2 (broad), so no change
        
        # Test expand_splits=True for serological to molecular
        a2_normal = self.converter.convert("A2", "m")
        a2_broad = self.converter.convert("A2", "m", expand_splits=True)
        
        # Should have more alleles with expand_splits
        self.assertGreater(len(a2_broad), len(a2_normal))
        
        # Normal A2 alleles should be included in expanded A2
        for allele in a2_normal:
            self.assertIn(allele, a2_broad)
        
        # A203 alleles should be included in expanded A2
        a203_alleles = self.converter.convert("A203", "m")
        for allele in a203_alleles:
            self.assertIn(allele, a2_broad)
        
        # A210 alleles should be included in expanded A2
        a210_alleles = self.converter.convert("A210", "m")
        for allele in a210_alleles:
            self.assertIn(allele, a2_broad)
    
    def test_broad_split_auto_detect(self):
        """Test broad/split functionality with auto-detection."""
        # Auto-detect with return_broad
        result = self.converter.convert("A*02:03", return_broad=True)
        self.assertEqual(result, "A2")
        
        # Auto-detect with expand_splits
        a2_normal = self.converter.convert("A2")
        a2_broad = self.converter.convert("A2", expand_splits=True)
        self.assertGreater(len(a2_broad), len(a2_normal))
    
    def test_broad_split_edge_cases(self):
        """Test edge cases for broad/split functionality."""
        # Test that non-split antigens aren't affected by return_broad
        result_normal = self.converter.convert("A*01:01", "s")
        result_broad = self.converter.convert("A*01:01", "s", return_broad=True)
        self.assertEqual(result_normal, result_broad)
        
        # Test that non-broad antigens aren't affected by expand_splits
        a1_normal = self.converter.convert("A1", "m")
        a1_broad = self.converter.convert("A1", "m", expand_splits=True)
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