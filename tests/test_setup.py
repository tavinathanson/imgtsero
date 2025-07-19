"""Test setup utilities for downloading real IMGT/HLA data."""

import os
import tempfile
import shutil
from imgtsero.downloader import download_data


class TestDataManager:
    """Manages test data download and cleanup."""
    
    def __init__(self):
        self.test_data_dir = None
        self.downloaded = False
    
    def setup_test_data(self, version="3610"):
        """Download real IMGT/HLA data for testing."""
        if self.downloaded:
            return self.test_data_dir
        
        self.test_data_dir = tempfile.mkdtemp(prefix="imgtsero_test_")
        
        try:
            print(f"Downloading real IMGT/HLA data version {version} for testing...")
            download_data(version, self.test_data_dir, verbose=False)
            self.downloaded = True
            print("Test data downloaded successfully.")
            return self.test_data_dir
        except Exception as e:
            print(f"Failed to download test data: {e}")
            self.cleanup()
            raise
    
    def cleanup(self):
        """Clean up test data directory."""
        if self.test_data_dir and os.path.exists(self.test_data_dir):
            shutil.rmtree(self.test_data_dir)
            self.test_data_dir = None
            self.downloaded = False


# Global test data manager
test_data_manager = TestDataManager()