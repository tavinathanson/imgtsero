"""Tests for download functionality."""

import unittest
import os
import tempfile
import shutil
from unittest.mock import patch
from imgtsero.downloader import download_file, download_data, extract_zip


class TestDownloader(unittest.TestCase):
    """Test cases for downloader."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
    
    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.temp_dir)
    
    def test_download_file_with_invalid_url(self):
        """Test download file with invalid URL."""
        filepath = os.path.join(self.temp_dir, "test_file.txt")
        # This should fail gracefully
        result = download_file("http://invalid-domain-that-does-not-exist-9999.com/test.txt", filepath)
        self.assertFalse(result)
    
    def test_download_data_creates_correct_paths(self):
        """Test that download_data creates correct file paths."""
        # Mock the actual download to test path logic
        with patch('imgtsero.downloader.download_file') as mock_download:
            mock_download.return_value = True
            
            result = download_data("3610", self.temp_dir, verbose=False)
            
            # Check that correct files are referenced
            self.assertIn('rel_dna_ser_file', result)
            self.assertIn('rel_ser_ser_file', result)
            
            # Check file paths are correct
            self.assertTrue(result['rel_dna_ser_file'].endswith('rel_dna_ser.3610.txt'))
            self.assertTrue(result['rel_ser_ser_file'].endswith('rel_ser_ser.3610.txt'))
            
            # Check that download was called with correct URLs
            calls = mock_download.call_args_list
            self.assertEqual(len(calls), 2)
            
            # Check URLs contain correct version and file names
            url1 = calls[0][0][0]
            url2 = calls[1][0][0]
            self.assertIn('3610', url1)
            self.assertIn('3610', url2)
            self.assertIn('rel_dna_ser.txt', url1)
            self.assertIn('rel_ser_ser.txt', url2)
    
    def test_extract_zip_nonexistent_file(self):
        """Test extract_zip with non-existent file."""
        fake_zip = os.path.join(self.temp_dir, "nonexistent.zip")
        result = extract_zip(fake_zip, self.temp_dir)
        self.assertFalse(result)
    
    def test_download_data_version_string_format(self):
        """Test download_data handles version string correctly."""
        with patch('imgtsero.downloader.download_file') as mock_download:
            mock_download.return_value = True
            
            # Test with different version formats
            result = download_data("3610", self.temp_dir, verbose=False)
            self.assertIsInstance(result, dict)
            
            # Version should be included in the URLs
            url_calls = [call[0][0] for call in mock_download.call_args_list]
            for url in url_calls:
                self.assertIn('/3610/', url)
    
    @patch('imgtsero.downloader.download_file')
    def test_download_data_failure_handling(self, mock_download_file):
        """Test download_data handles failures correctly."""
        # First file succeeds, second fails
        mock_download_file.side_effect = [True, False]
        
        with self.assertRaises(RuntimeError) as cm:
            download_data("3610", self.temp_dir, verbose=False)
        
        self.assertIn("Failed to download", str(cm.exception))


if __name__ == '__main__':
    unittest.main()