#!/usr/bin/env python
"""Unit test specifically for B*27:05 lookup fix."""

import unittest
import tempfile
import json
import os
from unittest.mock import patch, MagicMock
from imgtsero.kir_ligand import KIRLigandClassifier


class TestB27LookupFix(unittest.TestCase):
    """Test that B*27:05 lookup works after compression fix."""
    
    def setUp(self):
        """Set up test environment."""
        self.temp_dir = tempfile.mkdtemp()
        self.classifier = KIRLigandClassifier("3.61.0", self.temp_dir)
    
    def tearDown(self):
        """Clean up test environment."""
        import shutil
        shutil.rmtree(self.temp_dir, ignore_errors=True)
    
    def test_b27_05_found_after_compression(self):
        """Test that B*27:05 is found when only high-res alleles exist in API data."""
        # Simulate data with only high-resolution alleles
        test_data = {
            "B*27:05:02": "Bw4",
            "B*27:05:09": "Bw4", 
            "B*27:05:18": "Bw4",
            "B*27:05:02:01": "Bw4",
        }
        
        # Apply compression
        compressed = self.classifier._compress_to_four_digit(test_data)
        
        # B*27:05 should now be in the data
        self.assertIn("B*27:05", compressed)
        self.assertEqual(compressed["B*27:05"], "Bw4")
        
        # Save and load
        cache_file = os.path.join(self.temp_dir, "kir_ligand_3.61.0.json")
        with open(cache_file, 'w') as f:
            json.dump(compressed, f)
        
        self.classifier.load_data()
        
        # Test lookup works
        result = self.classifier.get_kir_ligand("B*27:05")
        self.assertEqual(result, "Bw4")
        
        # Test classification
        classification = self.classifier.classify_allele("B*27:05")
        self.assertEqual(classification["kir_ligand_type"], "Bw4")
        self.assertTrue(classification["is_kir_ligand"])
        self.assertEqual(classification["kir_receptors"], ["KIR3DL1"])
    
    @patch('urllib.request.urlopen')
    def test_real_world_api_scenario(self, mock_urlopen):
        """Test the real-world scenario where API returns only high-res B*27:05 alleles."""
        # Mock API response with realistic data
        api_response = {
            "data": [
                # B*27:05 alleles - only high resolution
                {"name": "B*27:05:02", "matching.kir_ligand": "Bw4"},
                {"name": "B*27:05:02:01", "matching.kir_ligand": "Bw4"},
                {"name": "B*27:05:09", "matching.kir_ligand": "Bw4"},
                {"name": "B*27:05:18", "matching.kir_ligand": "Bw4"},
                {"name": "B*27:05:26", "matching.kir_ligand": "Bw4"},
                
                # B*27:02 alleles - including 4-digit
                {"name": "B*27:02", "matching.kir_ligand": "Bw4"},
                {"name": "B*27:02:01", "matching.kir_ligand": "Bw4"},
                
                # Some non-Bw4 alleles
                {"name": "B*07:02", "matching.kir_ligand": "Bw6"},
                {"name": "B*07:02:01", "matching.kir_ligand": "Bw6"},
            ]
        }
        
        mock_response = MagicMock()
        mock_response.read.return_value = json.dumps(api_response).encode('utf-8')
        mock_urlopen.return_value.__enter__.return_value = mock_response
        
        # Force download from "API"
        self.classifier.load_data(force_download=True)
        
        # B*27:05 should work now
        self.assertEqual(self.classifier.get_kir_ligand("B*27:05"), "Bw4")
        
        # B*27:02 should still work (was already in data)
        self.assertEqual(self.classifier.get_kir_ligand("B*27:02"), "Bw4")
        
        # High-res lookups should still work
        self.assertEqual(self.classifier.get_kir_ligand("B*27:05:02"), "Bw4")
        self.assertEqual(self.classifier.get_kir_ligand("B*27:05:09"), "Bw4")
    
    def test_multiple_4digit_alleles_compressed(self):
        """Test that multiple different 4-digit alleles are compressed correctly."""
        test_data = {
            # B*27:05 family
            "B*27:05:02": "Bw4",
            "B*27:05:09": "Bw4",
            
            # B*44:02 family  
            "B*44:02:01": "Bw4",
            "B*44:02:01:01": "Bw4",
            
            # B*51:01 family
            "B*51:01:01": "Bw4",
            "B*51:01:01:01": "Bw4",
            
            # C*01:02 family
            "C*01:02:01": "C1",
            "C*01:02:29": "C1",
        }
        
        compressed = self.classifier._compress_to_four_digit(test_data)
        
        # Check all 4-digit forms were added
        self.assertIn("B*27:05", compressed)
        self.assertEqual(compressed["B*27:05"], "Bw4")
        
        self.assertIn("B*44:02", compressed)
        self.assertEqual(compressed["B*44:02"], "Bw4")
        
        self.assertIn("B*51:01", compressed)
        self.assertEqual(compressed["B*51:01"], "Bw4")
        
        self.assertIn("C*01:02", compressed)
        self.assertEqual(compressed["C*01:02"], "C1")


if __name__ == '__main__':
    unittest.main()