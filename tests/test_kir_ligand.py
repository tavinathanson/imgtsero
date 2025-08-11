"""Tests for KIR ligand classification functionality."""

import os
import json
import tempfile
import unittest
from unittest.mock import patch, MagicMock
import urllib.error

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
        
        # Verify results
        self.assertEqual(len(result), len(self.sample_data))
        for allele, expected_kir in self.sample_data.items():
            self.assertEqual(result.get(allele), expected_kir)
    
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
        self.assertEqual(result["kir_receptors"], ["KIR3DL1"])
        self.assertEqual(result["source"], "api")
        
        # Test Bw6 allele (not KIR ligand)
        result = self.classifier.classify_allele("B*07:02:01")
        self.assertFalse(result["is_kir_ligand"])
        self.assertEqual(result["kir_ligand_type"], "Bw6")
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


if __name__ == '__main__':
    unittest.main()