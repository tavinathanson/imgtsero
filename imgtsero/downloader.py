"""Download functionality for IMGT/HLA data files."""

import os
import urllib.request
import urllib.error
import zipfile


def download_file(url, filepath):
    """Download a file from URL to filepath."""
    try:
        urllib.request.urlretrieve(url, filepath)
        return True
    except urllib.error.URLError as e:
        print(f"Error downloading {url}: {e}")
        return False


def extract_zip(zip_path, extract_to):
    """Extract a ZIP file to the specified directory."""
    try:
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(extract_to)
        return True
    except zipfile.BadZipFile:
        print(f"Error: {zip_path} is not a valid ZIP file")
        return False
    except Exception as e:
        print(f"Error extracting {zip_path}: {e}")
        return False


def download_data(version, data_dir="data", verbose=True):
    """Download IMGT/HLA WMDA data files for the specified version."""
    # Create data directory
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
    
    # URLs for WMDA files
    base_url = f"https://raw.githubusercontent.com/ANHIG/IMGTHLA/{version}/wmda"
    rel_dna_ser_url = f"{base_url}/rel_dna_ser.txt"
    rel_ser_ser_url = f"{base_url}/rel_ser_ser.txt"
    
    # File paths
    rel_dna_ser_file = os.path.join(data_dir, f"rel_dna_ser.{version}.txt")
    rel_ser_ser_file = os.path.join(data_dir, f"rel_ser_ser.{version}.txt")
    
    # Download molecular to serological mapping
    if verbose:
        print(f"Downloading rel_dna_ser.{version}.txt...")
    if not download_file(rel_dna_ser_url, rel_dna_ser_file):
        raise RuntimeError(f"Failed to download rel_dna_ser.{version}.txt")
    if verbose:
        print(f"Downloaded rel_dna_ser.{version}.txt")
        print(f"Saved to: {rel_dna_ser_file}")
    
    # Download serological relationships mapping
    if verbose:
        print(f"Downloading rel_ser_ser.{version}.txt...")
    if not download_file(rel_ser_ser_url, rel_ser_ser_file):
        raise RuntimeError(f"Failed to download rel_ser_ser.{version}.txt")
    if verbose:
        print(f"Downloaded rel_ser_ser.{version}.txt")
        print(f"Saved to: {rel_ser_ser_file}")
        print("Download completed successfully!")
    
    return {
        'rel_dna_ser_file': rel_dna_ser_file,
        'rel_ser_ser_file': rel_ser_ser_file
    }