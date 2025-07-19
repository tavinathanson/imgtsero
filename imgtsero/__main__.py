#!/usr/bin/env python3
"""
IMGT/HLA data downloader

Downloads allele list and HLA data files from the official IMGT/HLA GitHub repository.
Usage: python -m imgtsero <version>
"""

import sys
from .downloader import download_data


def main():
    """Main function to download IMGT/HLA data files."""
    if len(sys.argv) != 2:
        print("Usage: python -m imgtsero <version>")
        print("Example: python -m imgtsero 3610")
        sys.exit(1)
    
    version = sys.argv[1]
    
    try:
        download_data(version)
    except RuntimeError as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()