"""Tests for KIR ligand classification functionality."""

import os
import json
import tempfile
import unittest
from unittest.mock import patch, MagicMock
import urllib.error
import pytest

from imgtsero.kir_ligand import KIRLigandClassifier
from imgtsero import HLAConverter


class TestKIRLigandClassifier(unittest.TestCase):
    """Test cases for KIRLigandClassifier."""
    
    def setUp(self):
        """Set up test environment."""
        self.temp_dir = tempfile.mkdtemp()
        self.test_version = "3.61.0"
        self.classifier = KIRLigandClassifier(self.test_version, self.temp_dir)
        
        # Sample KIR ligand data for testing
        self.sample_data = {
            "B*07:02:01": "Bw6",
            "B*07:02": "Bw6",
            "B*27:05:02": "Bw4",
            "B*27:05": "Bw4",
            "B*44:02:01:01": "Bw4",
            "B*44:02": "Bw4",
            "C*01:02:01": "C1",
            "C*01:02": "C1",
            "C*07:02:01:01": "C2",
            "C*07:02": "C2",
            "A*23:01:01": "Bw4",
            "A*23:01": "Bw4",
            "A*24:02:01:01": "Bw4",
            "A*24:02": "Bw4",
            "A*01:01:01:01": None,
            "A*01:01": None
        }
    
    def tearDown(self):
        """Clean up test environment."""
        # Remove temp directory and its contents
        import shutil
        shutil.rmtree(self.temp_dir, ignore_errors=True)
    
    def test_init(self):
        """Test classifier initialization."""
        self.assertEqual(self.classifier.version, self.test_version)
        self.assertEqual(self.classifier.data_dir, self.temp_dir)
        self.assertFalse(self.classifier._loaded)
    
    def test_version_normalization(self):
        """Test version format conversion."""
        # Test 3610 -> 3.61.0
        classifier1 = KIRLigandClassifier(3610, self.temp_dir)
        self.assertEqual(classifier1.version, "3.61.0")
        
        # Test "3610" -> 3.61.0
        classifier2 = KIRLigandClassifier("3610", self.temp_dir)
        self.assertEqual(classifier2.version, "3.61.0")
        
        # Test 3.61.0 stays same
        classifier3 = KIRLigandClassifier("3.61.0", self.temp_dir)
        self.assertEqual(classifier3.version, "3.61.0")
        
        # Test 3.61 -> 3.61.0
        classifier4 = KIRLigandClassifier("3.61", self.temp_dir)
        self.assertEqual(classifier4.version, "3.61.0")
    
    def test_kir_receptor_mappings(self):
        """Test KIR receptor mappings."""
        self.assertEqual(self.classifier.kir_receptors["Bw4"], ["KIR3DL1"])
        self.assertEqual(self.classifier.kir_receptors["C1"], ["KIR2DL2", "KIR2DL3"])
        self.assertEqual(self.classifier.kir_receptors["C2"], ["KIR2DL1"])
    
    @patch('urllib.request.urlopen')
    def test_fetch_kir_data_success(self, mock_urlopen):
        """Test successful API data fetch."""
        # Mock API response
        api_response = {
            "data": [
                {"name": allele, "matching.kir_ligand": kir_ligand}
                for allele, kir_ligand in self.sample_data.items()
            ]
        }
        
        mock_response = MagicMock()
        mock_response.read.return_value = json.dumps(api_response).encode('utf-8')
        mock_urlopen.return_value.__enter__.return_value = mock_response
        
        # Fetch data
        result = self.classifier._fetch_kir_data_from_api()
        
        # Verify results - after compression, we'll have more entries
        # Check that all original data is present
        for allele, expected_kir in self.sample_data.items():
            self.assertEqual(result.get(allele), expected_kir)
        
        # Check that 4-digit forms were added for consistent alleles
        # B*27:05 should be added from B*27:05:02
        self.assertIn("B*27:05", result)
        self.assertEqual(result["B*27:05"], "Bw4")
        
        # A*23 and A*24 should be added (2-digit A forms)
        self.assertIn("A*23", result)
        self.assertEqual(result["A*23"], "Bw4")
        self.assertIn("A*24", result)
        self.assertEqual(result["A*24"], "Bw4")
    
    @patch('urllib.request.urlopen')
    def test_fetch_kir_data_no_results(self, mock_urlopen):
        """Test API fetch with no results."""
        # Mock empty API response
        api_response = {"data": []}
        
        mock_response = MagicMock()
        mock_response.read.return_value = json.dumps(api_response).encode('utf-8')
        mock_urlopen.return_value.__enter__.return_value = mock_response
        
        # Should raise error for no data
        with self.assertRaises(RuntimeError) as cm:
            self.classifier._fetch_kir_data_from_api()
        self.assertIn("No KIR ligand data found", str(cm.exception))
    
    @patch('urllib.request.urlopen')
    def test_fetch_kir_data_http_error(self, mock_urlopen):
        """Test API fetch with HTTP error."""
        mock_urlopen.side_effect = urllib.error.HTTPError(None, 404, "Not Found", None, None)
        
        with self.assertRaises(RuntimeError) as cm:
            self.classifier._fetch_kir_data_from_api()
        self.assertIn("not available for version", str(cm.exception))
    
    def test_cache_operations(self):
        """Test cache save and load operations."""
        # Set test data
        self.classifier._kir_ligand_map = self.sample_data
        
        # Save to cache
        self.classifier._save_to_cache()
        
        # Verify cache file exists
        cache_file = self.classifier._get_cache_filename()
        self.assertTrue(os.path.exists(cache_file))
        
        # Create new classifier and load from cache
        new_classifier = KIRLigandClassifier(self.test_version, self.temp_dir)
        self.assertTrue(new_classifier._load_from_cache())
        
        # Verify loaded data matches
        self.assertEqual(new_classifier._kir_ligand_map, self.sample_data)
    
    def test_get_kir_ligand(self):
        """Test KIR ligand lookup."""
        # Set test data
        self.classifier._kir_ligand_map = self.sample_data
        self.classifier._loaded = True
        
        # Test exact match
        self.assertEqual(self.classifier.get_kir_ligand("B*07:02:01"), "Bw6")
        self.assertEqual(self.classifier.get_kir_ligand("B*27:05:02"), "Bw4")
        self.assertEqual(self.classifier.get_kir_ligand("C*01:02:01"), "C1")
        
        # Test shorter form lookup
        self.assertEqual(self.classifier.get_kir_ligand("B*07:02:01:02"), "Bw6")
        self.assertEqual(self.classifier.get_kir_ligand("B*27:05:02:03"), "Bw4")
        
        # Test no match
        self.assertIsNone(self.classifier.get_kir_ligand("B*99:99:99"))
        
        # Test null KIR ligand
        self.assertIsNone(self.classifier.get_kir_ligand("A*01:01:01:01"))
    
    def test_classify_allele(self):
        """Test allele classification."""
        # Set test data
        self.classifier._kir_ligand_map = self.sample_data
        self.classifier._loaded = True
        
        # Test Bw4 allele (KIR ligand)
        result = self.classifier.classify_allele("B*27:05:02")
        self.assertTrue(result["is_kir_ligand"])
        self.assertEqual(result["kir_ligand_type"], "Bw4")
        self.assertEqual(result["kir_ligand_type_detail"], "Bw4")  # From test data
        self.assertEqual(result["kir_receptors"], ["KIR3DL1"])
        self.assertEqual(result["source"], "api")
        
        # Test Bw6 allele (not KIR ligand)
        result = self.classifier.classify_allele("B*07:02:01")
        self.assertFalse(result["is_kir_ligand"])
        self.assertEqual(result["kir_ligand_type"], "Bw6")
        self.assertEqual(result["kir_ligand_type_detail"], "Bw6")  # From test data
        self.assertEqual(result["kir_receptors"], [])
        self.assertEqual(result["source"], "api")
        
        # Test C1 allele
        result = self.classifier.classify_allele("C*01:02:01")
        self.assertTrue(result["is_kir_ligand"])
        self.assertEqual(result["kir_ligand_type"], "C1")
        self.assertEqual(result["kir_receptors"], ["KIR2DL2", "KIR2DL3"])
        
        # Test with matching bead annotation (no error)
        result = self.classifier.classify_allele("B*07:02:01", bead_annotation="Bw6")
        self.assertFalse(result["is_kir_ligand"])
        self.assertEqual(result["kir_ligand_type"], "Bw6")
        self.assertEqual(result["source"], "api")  # Always API
        
        # Test with conflicting bead annotation (should raise error)
        with self.assertRaises(ValueError) as cm:
            self.classifier.classify_allele("B*27:05:02", bead_annotation="Bw6")
        self.assertIn("conflicts with API data", str(cm.exception))
    
    def test_classify_serological(self):
        """Test serological antigen classification."""
        # Set test data
        self.classifier._kir_ligand_map = self.sample_data
        self.classifier._loaded = True
        
        # Test consistent KIR ligand type
        molecular_alleles = ["B*07:02:01", "B*07:02:01:02"]
        result = self.classifier.classify_serological("B7", molecular_alleles)
        self.assertFalse(result["is_kir_ligand"])
        self.assertEqual(result["kir_ligand_type"], "Bw6")
        self.assertEqual(result["source"], "api")
        
        # Test with matching bead annotation
        result = self.classifier.classify_serological("B7", molecular_alleles, bead_annotation="Bw6")
        self.assertFalse(result["is_kir_ligand"])
        self.assertEqual(result["kir_ligand_type"], "Bw6")
        self.assertEqual(result["source"], "api")
        
        # Test with conflicting bead annotation
        with self.assertRaises(ValueError) as cm:
            self.classifier.classify_serological("B7", molecular_alleles, bead_annotation="Bw4")
        self.assertIn("conflicts with API data", str(cm.exception))
        
        # Test mixed types (shouldn't happen in practice)
        mixed_alleles = ["B*07:02:01", "B*27:05:02"]
        result = self.classifier.classify_serological("Mixed", mixed_alleles)
        self.assertFalse(result["is_kir_ligand"])
        self.assertIn("Mixed:", result["kir_ligand_type"])
        
        # Test no KIR ligand data
        no_kir_alleles = ["B*99:99:99"]
        result = self.classifier.classify_serological("Unknown", no_kir_alleles)
        self.assertFalse(result["is_kir_ligand"])
        self.assertIsNone(result["kir_ligand_type"])
    
    def test_get_all_kir_ligands(self):
        """Test getting all KIR ligands grouped by type."""
        # Set test data
        self.classifier._kir_ligand_map = self.sample_data
        self.classifier._loaded = True
        
        grouped = self.classifier.get_all_kir_ligands()
        
        # Check Bw4 group
        self.assertIn("Bw4", grouped)
        bw4_alleles = grouped["Bw4"]
        self.assertIn("B*27:05:02", bw4_alleles)
        self.assertIn("B*44:02:01:01", bw4_alleles)
        self.assertIn("A*23:01:01", bw4_alleles)
        
        # Check Bw6 group
        self.assertIn("Bw6", grouped)
        bw6_alleles = grouped["Bw6"]
        self.assertIn("B*07:02:01", bw6_alleles)
        
        # Check C1 and C2 groups
        self.assertIn("C1", grouped)
        self.assertIn("C*01:02:01", grouped["C1"])
        self.assertIn("C2", grouped)
        self.assertIn("C*07:02:01:01", grouped["C2"])


class TestHLAConverterKIRIntegration(unittest.TestCase):
    """Test KIR ligand integration with HLAConverter."""
    
    def setUp(self):
        """Set up test environment."""
        self.temp_dir = tempfile.mkdtemp()
        
        # Create minimal test data files for HLAParser
        rel_dna_ser_file = os.path.join(self.temp_dir, "rel_dna_ser.3610.txt")
        with open(rel_dna_ser_file, 'w') as f:
            f.write("# Test data\n")
            f.write("B*;07:02:01;7;;;\n")
            f.write("B*;27:05:02;27;;;\n")
            f.write("C*;01:02:01;1;;;\n")
            f.write("C*;07:02:01:01;7;;;\n")
        
        rel_ser_ser_file = os.path.join(self.temp_dir, "rel_ser_ser.3610.txt")
        with open(rel_ser_ser_file, 'w') as f:
            f.write("# Test data\n")
    
    def tearDown(self):
        """Clean up test environment."""
        import shutil
        shutil.rmtree(self.temp_dir, ignore_errors=True)
    
    def test_converter_without_kir(self):
        """Test converter without KIR enabled."""
        converter = HLAConverter(version=3610, data_dir=self.temp_dir, enable_kir=False)
        self.assertIsNone(converter.kir_classifier)
        
        # Should raise error when trying to classify
        with self.assertRaises(RuntimeError) as cm:
            converter.classify_kir_ligand("B*07:02:01")
        self.assertIn("KIR ligand classification not enabled", str(cm.exception))
    
    @patch.object(KIRLigandClassifier, 'load_data')
    @patch.object(KIRLigandClassifier, 'classify_allele')
    def test_converter_with_kir(self, mock_classify, mock_load):
        """Test converter with KIR enabled."""
        # Set up mocks
        mock_classify.return_value = {
            "is_kir_ligand": False,
            "kir_ligand_type": "Bw6",
            "kir_receptors": [],
            "source": "api"
        }
        
        converter = HLAConverter(version=3610, data_dir=self.temp_dir, enable_kir=True)
        self.assertIsNotNone(converter.kir_classifier)
        
        # Test molecular allele classification
        result = converter.classify_kir_ligand("B*07:02:01")
        mock_load.assert_called_once()
        mock_classify.assert_called_once_with("B*07:02:01", None)
        self.assertEqual(result["kir_ligand_type"], "Bw6")
        
        # Test with bead annotation
        mock_classify.reset_mock()
        result = converter.classify_kir_ligand("B*07:02:01", bead_annotation="Bw6")
        mock_classify.assert_called_once_with("B*07:02:01", "Bw6")


class TestKIRLigandRealFunctionality(unittest.TestCase):
    """Test KIR ligand functionality with real data (no mocks)."""
    
    def setUp(self):
        """Set up test environment."""
        self.temp_dir = tempfile.mkdtemp()
        self.classifier = KIRLigandClassifier("3.61.0", self.temp_dir)
        
        # Create a test cache file with known data
        test_data = {
            "B*07:02:01": "Bw6",
            "B*07:02": "Bw6",
            "B*27:05:02": "Bw4",
            "B*27:05": "Bw4",
            "B*44:02:01:01": "Bw4",
            "B*44:02": "Bw4",
            "C*01:02:01": "C1",
            "C*01:02": "C1",
            "C*07:02:01:01": "C2",
            "C*07:02": "C2",
            "C*03:04:01": "C1",
            "C*03:04": "C1",
            "C*04:01:01": "C2",
            "C*04:01": "C2",
            "A*23:01:01": "Bw4",
            "A*23:01": "Bw4",
            "A*24:02:01:01": "Bw4",
            "A*24:02": "Bw4",
            "A*32:01:01": "Bw4",
            "A*01:01:01:01": None,
            "A*01:01": None,
            "A*02:01:01:01": None,
            "B*08:01:01": "Bw6",
            "B*15:01:01:01": "Bw6",
            "B*35:01:01": "Bw6",
            "B*13:02:01": "Bw4",
            "B*37:01:01": "Bw4",
            "B*38:01:01": "Bw6",
            "B*39:01:01": "Bw6",
            "B*51:01:01": "Bw4",
            "B*52:01:01": "Bw4",
            "B*53:01:01": "Bw4",
            "B*57:01:01": "Bw4",
            "B*58:01:01": "Bw4"
        }
        
        # Save test data to cache
        cache_file = os.path.join(self.temp_dir, "kir_ligand_3.61.0.json")
        with open(cache_file, 'w') as f:
            json.dump(test_data, f)
    
    def tearDown(self):
        """Clean up test environment."""
        import shutil
        shutil.rmtree(self.temp_dir, ignore_errors=True)
    
    def test_real_classification_flow(self):
        """Test complete classification flow with cached data."""
        # Load from cache (not API)
        self.classifier.load_data()
        
        # Test various B alleles
        b_tests = [
            ("B*07:02:01", "Bw6", False),
            ("B*27:05:02", "Bw4", True),
            ("B*44:02:01:01", "Bw4", True),
            ("B*08:01:01", "Bw6", False),
            ("B*51:01:01", "Bw4", True),
            ("B*57:01:01", "Bw4", True),
        ]
        
        for allele, expected_type, expected_ligand in b_tests:
            result = self.classifier.classify_allele(allele)
            self.assertEqual(result["kir_ligand_type"], expected_type)
            self.assertEqual(result["is_kir_ligand"], expected_ligand)
            if expected_ligand and expected_type == "Bw4":
                self.assertEqual(result["kir_receptors"], ["KIR3DL1"])
    
    def test_real_c_allele_classification(self):
        """Test C allele classification with real data."""
        self.classifier.load_data()
        
        # Test C1 alleles
        c1_alleles = ["C*01:02:01", "C*03:04:01"]
        for allele in c1_alleles:
            result = self.classifier.classify_allele(allele)
            self.assertEqual(result["kir_ligand_type"], "C1")
            self.assertTrue(result["is_kir_ligand"])
            self.assertEqual(result["kir_receptors"], ["KIR2DL2", "KIR2DL3"])
        
        # Test C2 alleles
        c2_alleles = ["C*07:02:01:01", "C*04:01:01"]
        for allele in c2_alleles:
            result = self.classifier.classify_allele(allele)
            self.assertEqual(result["kir_ligand_type"], "C2")
            self.assertTrue(result["is_kir_ligand"])
            self.assertEqual(result["kir_receptors"], ["KIR2DL1"])
    
    def test_real_a_allele_classification(self):
        """Test A allele classification with real data."""
        self.classifier.load_data()
        
        # Test Bw4-positive A alleles
        bw4_a_alleles = ["A*23:01:01", "A*24:02:01:01"]
        for allele in bw4_a_alleles:
            result = self.classifier.classify_allele(allele)
            self.assertEqual(result["kir_ligand_type"], "Bw4")
            self.assertTrue(result["is_kir_ligand"])
            self.assertEqual(result["kir_receptors"], ["KIR3DL1"])
        
        # Test non-KIR ligand A alleles
        non_kir_a_alleles = ["A*01:01:01:01", "A*02:01:01:01"]
        for allele in non_kir_a_alleles:
            result = self.classifier.classify_allele(allele)
            self.assertIsNone(result["kir_ligand_type"])
            self.assertFalse(result["is_kir_ligand"])
            self.assertEqual(result["kir_receptors"], [])
    
    def test_progressive_resolution(self):
        """Test that shorter allele forms are found when exact match fails."""
        self.classifier.load_data()
        
        # These don't exist in our test data but shorter forms do
        test_cases = [
            ("B*07:02:01:02:03", "Bw6"),  # Should find B*07:02:01 or B*07:02
            ("B*27:05:02:01", "Bw4"),      # Should find B*27:05:02 or B*27:05
            ("C*01:02:01:02", "C1"),        # Should find C*01:02:01 or C*01:02
        ]
        
        for allele, expected_type in test_cases:
            result = self.classifier.get_kir_ligand(allele)
            self.assertEqual(result, expected_type)
    
    def test_validation_with_bead_annotations(self):
        """Test bead annotation validation against API data."""
        self.classifier.load_data()
        
        # Test matching annotations (should not raise)
        result = self.classifier.classify_allele("B*07:02:01", bead_annotation="Bw6")
        self.assertEqual(result["kir_ligand_type"], "Bw6")
        
        result = self.classifier.classify_allele("B*27:05:02", bead_annotation="Bw4")
        self.assertEqual(result["kir_ligand_type"], "Bw4")
        
        # Test conflicting annotations (should raise)
        with self.assertRaises(ValueError) as cm:
            self.classifier.classify_allele("B*07:02:01", bead_annotation="Bw4")
        self.assertIn("conflicts with API data", str(cm.exception))
        
        with self.assertRaises(ValueError) as cm:
            self.classifier.classify_allele("B*27:05:02", bead_annotation="Bw6")
        self.assertIn("conflicts with API data", str(cm.exception))
    
    def test_parse_bead_annotation(self):
        """Test parsing of SAB bead names."""
        test_cases = [
            ("B27,Bw4", ("B27", "Bw4")),
            ("B7,Bw6", ("B7", "Bw6")),
            ("Cw7", ("Cw7", None)),
            ("C1", ("C1", None)),
            ("B44,Bw4", ("B44", "Bw4")),
            ("A23,Bw4", ("A23", "Bw4")),
            ("B27,something", ("B27", None)),  # Invalid annotation
            ("B27", ("B27", None)),            # No annotation
        ]
        
        for bead_name, expected in test_cases:
            result = KIRLigandClassifier.parse_bead_annotation(bead_name)
            self.assertEqual(result, expected)
    
    def test_integration_with_converter(self):
        """Test full integration with HLAConverter."""
        # Create test WMDA data
        rel_dna_ser_file = os.path.join(self.temp_dir, "rel_dna_ser.3610.txt")
        with open(rel_dna_ser_file, 'w') as f:
            f.write("# Test data\n")
            f.write("B*;07:02:01;7;;;\n")
            f.write("B*;07:02:01:01;7;;;\n")
            f.write("B*;27:05:02;27;;;\n")
            f.write("B*;27:05:02:01;27;;;\n")
            f.write("C*;01:02:01;1;;;\n")
            f.write("C*;07:02:01:01;7;;;\n")
            f.write("A*;23:01:01;23;;;\n")
            f.write("A*;01:01:01:01;1;;;\n")
        
        rel_ser_ser_file = os.path.join(self.temp_dir, "rel_ser_ser.3610.txt")
        with open(rel_ser_ser_file, 'w') as f:
            f.write("# Test data\n")
        
        # Initialize converter with KIR
        converter = HLAConverter(version=3610, data_dir=self.temp_dir, enable_kir=True)
        
        # Test molecular classification
        result = converter.classify_kir_ligand("B*07:02:01")
        self.assertEqual(result["kir_ligand_type"], "Bw6")
        self.assertFalse(result["is_kir_ligand"])
        
        # Test serological classification
        result = converter.classify_kir_ligand("B27")
        self.assertEqual(result["kir_ligand_type"], "Bw4")
        self.assertTrue(result["is_kir_ligand"])
        
        # Test with validation
        result = converter.classify_kir_ligand("B7", bead_annotation="Bw6")
        self.assertEqual(result["kir_ligand_type"], "Bw6")
        
        # Test validation error
        with self.assertRaises(ValueError):
            converter.classify_kir_ligand("B27", bead_annotation="Bw6")
    
    def test_cache_persistence(self):
        """Test that cache is properly saved and loaded."""
        # Load data (from our pre-created cache)
        self.classifier.load_data()
        
        # Create a new classifier instance
        new_classifier = KIRLigandClassifier("3.61.0", self.temp_dir)
        
        # Should load from cache without API call
        new_classifier.load_data()
        
        # Verify data is loaded correctly
        result = new_classifier.classify_allele("B*27:05:02")
        self.assertEqual(result["kir_ligand_type"], "Bw4")
        self.assertTrue(result["is_kir_ligand"])
    
    def test_grouped_kir_ligands(self):
        """Test grouping of alleles by KIR ligand type."""
        self.classifier.load_data()
        
        grouped = self.classifier.get_all_kir_ligands()
        
        # Verify groups exist
        self.assertIn("Bw4", grouped)
        self.assertIn("Bw6", grouped)
        self.assertIn("C1", grouped)
        self.assertIn("C2", grouped)
        
        # Verify some specific alleles are in correct groups
        self.assertIn("B*27:05:02", grouped["Bw4"])
        self.assertIn("B*44:02:01:01", grouped["Bw4"])
        self.assertIn("A*23:01:01", grouped["Bw4"])
        self.assertIn("B*07:02:01", grouped["Bw6"])
        self.assertIn("C*01:02:01", grouped["C1"])
        self.assertIn("C*07:02:01:01", grouped["C2"])
        
        # Verify counts are reasonable
        self.assertGreater(len(grouped["Bw4"]), 5)
        self.assertGreater(len(grouped["Bw6"]), 5)
        self.assertGreater(len(grouped["C1"]), 1)
        self.assertGreater(len(grouped["C2"]), 1)


class TestKIRLigandHighResolutionIssue(unittest.TestCase):
    """Test case to demonstrate the issue with B27*05 not returning Bw4."""
    
    def setUp(self):
        """Set up test environment."""
        self.temp_dir = tempfile.mkdtemp()
        self.classifier = KIRLigandClassifier("3.61.0", self.temp_dir)
        
        # Simulate API data that only has high-resolution alleles, not 4-digit ones
        # This mimics the real issue where B*27:05 is missing but B*27:05:09 exists
        test_data = {
            # High resolution B*27 alleles with Bw4
            "B*27:05:02": "Bw4",
            "B*27:05:09": "Bw4",
            "B*27:05:18": "Bw4",
            "B*27:05:02:01": "Bw4",
            # Note: B*27:05 is missing!
            
            # Other alleles for comparison
            "B*07:02": "Bw6",
            "B*07:02:01": "Bw6",
        }
        
        # Save test data to cache
        cache_file = os.path.join(self.temp_dir, "kir_ligand_3.61.0.json")
        with open(cache_file, 'w') as f:
            json.dump(test_data, f)
    
    def tearDown(self):
        """Clean up test environment."""
        import shutil
        shutil.rmtree(self.temp_dir, ignore_errors=True)
    
    def test_b27_05_lookup_fails_without_compression(self):
        """Test that B*27:05 lookup fails when only higher resolution data exists."""
        self.classifier.load_data()
        
        # This should fail because B*27:05 is not in the data
        result = self.classifier.get_kir_ligand("B*27:05")
        self.assertIsNone(result, "B*27:05 should not be found without compression")
        
        # But higher resolution lookups work
        result = self.classifier.get_kir_ligand("B*27:05:02")
        self.assertEqual(result, "Bw4")
        
        result = self.classifier.get_kir_ligand("B*27:05:09")
        self.assertEqual(result, "Bw4")
    
    def test_progressive_resolution_doesnt_help_for_missing_4digit(self):
        """Test that progressive resolution doesn't help when looking up 4-digit alleles."""
        self.classifier.load_data()
        
        # Looking up B*27:05:99 (doesn't exist) won't find B*27:05 because B*27:05 isn't in the data
        result = self.classifier.get_kir_ligand("B*27:05:99")
        self.assertIsNone(result, "Progressive resolution can't find B*27:05 if it's not in the data")


class TestKIRLigandCompression(unittest.TestCase):
    """Test the compression functionality for KIR ligand data."""
    
    def setUp(self):
        """Set up test environment."""
        self.temp_dir = tempfile.mkdtemp()
        self.classifier = KIRLigandClassifier("3.61.0", self.temp_dir)
    
    def tearDown(self):
        """Clean up test environment."""
        import shutil
        shutil.rmtree(self.temp_dir, ignore_errors=True)
    
    def test_compression_adds_four_digit_forms(self):
        """Test that compression adds 4-digit forms when consistent."""
        # Test data with only high-resolution alleles
        test_data = {
            "B*27:05:02": "Bw4",
            "B*27:05:09": "Bw4",
            "B*27:05:18": "Bw4",
            "B*07:02:01": "Bw6",
            "B*07:02:01:01": "Bw6",
            "C*01:02:01": "C1",
            "C*01:02:29": "C1",
        }
        
        # Apply compression
        compressed = self.classifier._compress_to_four_digit(test_data)
        
        # Check that 4-digit forms were added
        self.assertIn("B*27:05", compressed)
        self.assertEqual(compressed["B*27:05"], "Bw4")
        
        self.assertIn("B*07:02", compressed)
        self.assertEqual(compressed["B*07:02"], "Bw6")
        
        self.assertIn("C*01:02", compressed)
        self.assertEqual(compressed["C*01:02"], "C1")
        
        # Original alleles should still be present
        self.assertIn("B*27:05:02", compressed)
        self.assertIn("B*07:02:01", compressed)
    
    def test_compression_detects_inconsistencies(self):
        """Test that compression detects and reports inconsistencies."""
        # Test data with conflicting KIR ligand types for same 4-digit allele
        test_data = {
            "B*27:05:02": "Bw4",
            "B*27:05:09": "Bw6",  # Conflict!
            "B*07:02:01": "Bw6",
        }
        
        # Should raise ValueError for inconsistency
        with self.assertRaises(ValueError) as cm:
            self.classifier._compress_to_four_digit(test_data)
        
        error_msg = str(cm.exception)
        self.assertIn("Inconsistent KIR ligand types", error_msg)
        self.assertIn("B*27:05", error_msg)
        self.assertIn("B*27:05:02: Bw4", error_msg)
        self.assertIn("B*27:05:09: Bw6", error_msg)
    
    def test_compression_handles_null_values(self):
        """Test that compression correctly handles null KIR ligand values."""
        test_data = {
            "A*01:01:01:01": None,
            "A*01:01:01:02": None,
            "B*27:05:02": "Bw4",
            "B*27:05:09": "Bw4",
        }
        
        compressed = self.classifier._compress_to_four_digit(test_data)
        
        # B*27:05 should be added
        self.assertIn("B*27:05", compressed)
        self.assertEqual(compressed["B*27:05"], "Bw4")
        
        # A*01:01 should NOT be added (null values are skipped)
        self.assertNotIn("A*01:01", compressed)
        
        # Original null entries should still be present
        self.assertIn("A*01:01:01:01", compressed)
        self.assertIsNone(compressed["A*01:01:01:01"])
    
    def test_compression_handles_a_locus_two_digit(self):
        """Test that compression adds 2-digit forms for A locus alleles."""
        test_data = {
            "A*23:01:01": "Bw4",
            "A*23:01:01:01": "Bw4",
            "A*24:02:01:01": "Bw4",
            "A*24:02:01:02": "Bw4",
        }
        
        compressed = self.classifier._compress_to_four_digit(test_data)
        
        # Check that 2-digit A forms were added
        self.assertIn("A*23", compressed)
        self.assertEqual(compressed["A*23"], "Bw4")
        
        self.assertIn("A*24", compressed)
        self.assertEqual(compressed["A*24"], "Bw4")
    
    def test_compression_skips_none_values(self):
        """Test that compression ignores None values when checking consistency."""
        test_data = {
            # C*17:01 family - one None, others C2
            "C*17:01:01:01": None,  # No KIR data
            "C*17:01:01:02": "C2",
            "C*17:01:02": "C2",
            "C*17:01:03": "C2",
            # B*44:02 family - all have data
            "B*44:02:01:01": "Bw4 - 80T",
            "B*44:02:01:02": "Bw4 - 80T",
            # C*99:01 family - all None
            "C*99:01:01": None,
            "C*99:01:02": None,
        }
        
        compressed = self.classifier._compress_to_four_digit(test_data)
        
        # C*17:01 should be compressed to C2 (ignoring the None)
        self.assertIn("C*17:01", compressed)
        self.assertEqual(compressed["C*17:01"], "C2")
        
        # B*44:02 should be compressed normally
        self.assertIn("B*44:02", compressed)
        self.assertEqual(compressed["B*44:02"], "Bw4 - 80T")
        
        # C*99:01 should NOT be compressed (all None)
        self.assertNotIn("C*99:01", compressed)
        
        # Original entries should still be present
        self.assertIn("C*17:01:01:01", compressed)
        self.assertIsNone(compressed["C*17:01:01:01"])
    
    def test_compression_mixed_none_inconsistent(self):
        """Test compression with mixed None values and inconsistent non-None values."""
        test_data = {
            # Mixed with inconsistency
            "B*98:01:01": None,
            "B*98:01:02": "Bw4",
            "B*98:01:03": "Bw6",  # Inconsistent!
        }
        
        # Should raise error for inconsistent non-None values
        with self.assertRaises(ValueError) as context:
            self.classifier._compress_to_four_digit(test_data)
        
        self.assertIn("B*98:01", str(context.exception))
        self.assertIn("Bw4", str(context.exception))
        self.assertIn("Bw6", str(context.exception))
    
    def test_c17_01_classification_with_real_data(self):
        """Test that C*17:01 gets classified correctly with real API data."""
        # This test uses the real cached data which should include full pagination
        classifier = KIRLigandClassifier(version="3.61.0")
        classifier.load_data()
        
        # Test C*17:01 - should be compressed to C2
        result = classifier.get_kir_ligand("C*17:01")
        self.assertEqual(result, "C2", "C*17:01 should be classified as C2 after compression")
        
        # Also verify the full classification
        converter = HLAConverter(version=3610, enable_kir=True)
        classification = converter.classify_kir_ligand("C*17:01")
        self.assertTrue(classification['is_kir_ligand'])
        self.assertEqual(classification['kir_ligand_type'], 'C2')
        self.assertEqual(classification['source'], 'api')
    
    def test_api_data_completeness(self):
        """Test that API returns expected number of alleles.
        
        Based on IPD-IMGT/HLA version 3.61.0 statistics:
        - Total HLA alleles: 43,225 (all loci)
        - Our API query fetches only A*, B*, C* loci for KIR ligand classification
        - Expected A+B+C total from API: 28,222
        
        Source: https://www.ebi.ac.uk/ipd/imgt/hla/about/statistics/growth/
        """
        classifier = KIRLigandClassifier(version="3.61.0")
        classifier.load_data()
        
        # Check total number of entries in cache
        total_entries = len(classifier._kir_ligand_map)
        
        # After compression, we expect more entries than raw API count
        # API returns 28,222 but compression adds 4-digit forms
        self.assertGreaterEqual(total_entries, 28222, 
            f"Expected at least 28,222 alleles from API, got {total_entries}")
        self.assertLessEqual(total_entries, 35000,
            f"Expected no more than 35,000 entries after compression, got {total_entries}")
        
        # Check specific counts for each locus
        a_count = sum(1 for k in classifier._kir_ligand_map.keys() if k.startswith('A*'))
        b_count = sum(1 for k in classifier._kir_ligand_map.keys() if k.startswith('B*'))
        c_count = sum(1 for k in classifier._kir_ligand_map.keys() if k.startswith('C*'))
        
        # These are the exact API counts for version 3.61.0:
        # A*: 8,746 (from API meta.total)
        # B*: 10,628 (from API meta.total)
        # C*: 8,848 (from API meta.total)
        # Total A+B+C: 28,222
        # After compression, we get additional 4-digit and 2-digit forms
        self.assertGreaterEqual(a_count, 8746, f"Expected at least 8,746 A alleles, got {a_count}")
        self.assertGreaterEqual(b_count, 10628, f"Expected at least 10,628 B alleles, got {b_count}")
        self.assertGreaterEqual(c_count, 8848, f"Expected at least 8,848 C alleles, got {c_count}")
        
        # Verify our pagination fetched the expected total
        # Raw A+B+C from API should be ~28,222 before compression
        raw_total = a_count + b_count + c_count
        self.assertGreaterEqual(raw_total, 28222,
            f"Expected at least 28,222 A+B+C alleles total, got {raw_total}")
        
        # Verify we have multiple C*17:01 variants (not just C*17:01:01:01)
        c17_01_variants = [k for k in classifier._kir_ligand_map.keys() 
                          if k.startswith('C*17:01') and len(k.split(':')) >= 3]
        self.assertGreaterEqual(len(c17_01_variants), 20, 
            f"Expected multiple C*17:01 variants, got {len(c17_01_variants)}")
    
    @pytest.mark.slow
    def test_total_database_alleles_by_pagination_v3_61_0(self):
        """Test total allele count by actually paginating through ALL alleles.
        
        This test fetches ALL HLA alleles (not just A/B/C) and verifies:
        1. The pagination count matches the API's meta.total field
        2. The total is close to the official statistics (43,225 as of July 2025)
        
        To run this test: pytest -m slow
        """
        import urllib.request
        import urllib.parse
        import time
        
        version = "3.61.0"
        base_url = "https://www.ebi.ac.uk/cgi-bin/ipd/api/allele"
        
        # Query for ALL alleles in this version (no locus filter)
        params = {
            'fields': 'name',
            'query': f'eq(release_version,"{version}")',
            'limit': '1000',
            'format': 'json'
        }
        
        all_alleles = set()
        next_url = None
        page_count = 0
        
        while True:
            if next_url:
                url = f"{base_url}{next_url}"
            else:
                query_string = urllib.parse.urlencode(params)
                url = f"{base_url}?{query_string}"
            
            with urllib.request.urlopen(url) as response:
                data = json.loads(response.read().decode('utf-8'))
            
            if data.get('data'):
                for entry in data['data']:
                    allele_name = entry.get('name', '')
                    if allele_name:
                        all_alleles.add(allele_name)
            
            page_count += 1
            next_url = data.get('meta', {}).get('next')
            if not next_url:
                break
            
            time.sleep(0.05)  # Be nice to the API
        
        actual_count = len(all_alleles)
        
        # Verify reasonable range (stats page shows 43,225)
        self.assertGreaterEqual(actual_count, 43000,
            f"Expected at least 43,000 total alleles, got {actual_count}")
        self.assertLessEqual(actual_count, 45000,
            f"Expected no more than 45,000 total alleles, got {actual_count}")
        
        # Verify against meta.total
        params['limit'] = '1'
        query_string = urllib.parse.urlencode(params)
        url = f"{base_url}?{query_string}"
        
        with urllib.request.urlopen(url) as response:
            data = json.loads(response.read().decode('utf-8'))
        
        meta_total = data.get('meta', {}).get('total', 0)
        
        # Pagination count should match meta.total exactly
        self.assertEqual(actual_count, meta_total,
            f"Pagination count ({actual_count}) should match meta.total ({meta_total})")
        
        # Verify we paginated through many pages
        self.assertGreaterEqual(page_count, 40,
            f"Expected at least 40 pages with limit=1000, got {page_count}")
    
    @patch('urllib.request.urlopen')
    def test_fetch_with_compression_integration(self, mock_urlopen):
        """Test that fetch API integrates compression correctly."""
        # Mock API response with only high-resolution alleles
        api_response = {
            "data": [
                {"name": "B*27:05:02", "matching.kir_ligand": "Bw4"},
                {"name": "B*27:05:09", "matching.kir_ligand": "Bw4"},
                {"name": "B*27:05:18", "matching.kir_ligand": "Bw4"},
                {"name": "B*07:02:01", "matching.kir_ligand": "Bw6"},
                {"name": "B*07:02:01:01", "matching.kir_ligand": "Bw6"},
            ]
        }
        
        mock_response = MagicMock()
        mock_response.read.return_value = json.dumps(api_response).encode('utf-8')
        mock_urlopen.return_value.__enter__.return_value = mock_response
        
        # Load data (which should fetch and compress)
        self.classifier.load_data(force_download=True)
        
        # Now B*27:05 should be found!
        result = self.classifier.get_kir_ligand("B*27:05")
        self.assertEqual(result, "Bw4", "B*27:05 should be found after compression")
        
        # And B*07:02 too
        result = self.classifier.get_kir_ligand("B*07:02")
        self.assertEqual(result, "Bw6", "B*07:02 should be found after compression")


if __name__ == '__main__':
    unittest.main()