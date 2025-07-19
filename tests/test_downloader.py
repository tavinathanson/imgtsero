"""Tests for download functionality."""

import unittest
import os
import tempfile
import shutil
from unittest.mock import patch, MagicMock
from imgtsero.downloader import download_file, download_data


class TestDownloader(unittest.TestCase):
    """Test cases for downloader."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
    
    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.temp_dir)
    
    @patch('imgtsero.downloader.urllib.request.urlretrieve')
    def test_download_file_success(self, mock_urlretrieve):
        """Test successful file download."""
        mock_urlretrieve.return_value = None
        
        filepath = os.path.join(self.temp_dir, "test_file.txt")
        result = download_file("http://example.com/test.txt", filepath)
        
        self.assertTrue(result)
        mock_urlretrieve.assert_called_once_with("http://example.com/test.txt", filepath)
    
    @patch('imgtsero.downloader.urllib.request.urlretrieve')
    def test_download_file_failure(self, mock_urlretrieve):
        """Test failed file download."""
        from urllib.error import URLError
        mock_urlretrieve.side_effect = URLError("Network error")
        
        filepath = os.path.join(self.temp_dir, "test_file.txt")
        result = download_file("http://example.com/test.txt", filepath)
        
        self.assertFalse(result)
    
    @patch('imgtsero.downloader.download_file')
    @patch('imgtsero.downloader.extract_zip')
    def test_download_data_success(self, mock_extract_zip, mock_download_file):
        """Test successful data download."""
        mock_download_file.return_value = True
        mock_extract_zip.return_value = True
        
        # Mock os.remove to avoid FileNotFoundError
        with patch('os.remove'):
            result = download_data("3610", self.temp_dir, verbose=False)
        
        self.assertIn('rel_dna_ser_file', result)
        self.assertIn('rel_ser_ser_file', result)
        # Downloading two WMDA files
        self.assertEqual(mock_download_file.call_count, 2)
        mock_extract_zip.assert_not_called()
    
    @patch('imgtsero.downloader.download_file')
    def test_download_data_failure(self, mock_download_file):
        """Test failed data download."""
        mock_download_file.return_value = False
        
        with self.assertRaises(RuntimeError):
            download_data("3610", self.temp_dir, verbose=False)


if __name__ == '__main__':
    unittest.main()