#!/usr/bin/env python3
"""Example usage of imgtsero library."""

import imgtsero

def demonstrate_conversion():
    """Demonstrate bidirectional HLA conversion functionality."""
    print("=== IMGTSERO Demo ===\n")
    
    # Test conversions with new simplified API
    test_cases = [
        ("A*01:01", None),      # Auto-detect: molecular -> serological
        ("A1", None),           # Auto-detect: serological -> molecular
        ("A*01:01", "s"),       # Explicit: molecular -> serological
        ("A1", "m"),            # Explicit: serological -> molecular
        ("B*27:05", "s"),       # Molecular -> serological
        ("B27", "m"),           # Serological -> molecular
    ]
    
    print("HLA Format Conversion Examples:")
    print("-" * 50)
    
    for hla_input, target_format in test_cases:
        try:
            result = imgtsero.convert(hla_input, target_format)
            format_desc = "auto-detect" if target_format is None else f"'{target_format}'"
            print(f"{hla_input:12} -> {format_desc:12} : {result}")
        except Exception as e:
            format_desc = "auto-detect" if target_format is None else f"'{target_format}'"
            print(f"{hla_input:12} -> {format_desc:12} : Error - {e}")
    
    print("\n=== Error Handling ===")
    error_cases = [
        ("A*99:99", None),      # Invalid molecular allele
        ("A999", None),         # Invalid serological allele
        ("not-hla", None),      # Invalid format
    ]
    
    for hla_input, target_format in error_cases:
        try:
            result = imgtsero.convert(hla_input, target_format)
            format_desc = "auto-detect" if target_format is None else f"'{target_format}'"
            print(f"{hla_input:12} -> {format_desc:12} : {result}")
        except imgtsero.HLAConversionError as e:
            format_desc = "auto-detect" if target_format is None else f"'{target_format}'"
            print(f"{hla_input:12} -> {format_desc:12} : ✗ {e}")
        except Exception as e:
            format_desc = "auto-detect" if target_format is None else f"'{target_format}'"
            print(f"{hla_input:12} -> {format_desc:12} : Error - {e}")
    
    print("\n=== Download Example ===")
    print("To download IMGT/HLA data:")
    print("  python -m imgtsero 3610")
    print("or")
    print("  imgtsero.download_data(3610)")


if __name__ == "__main__":
    demonstrate_conversion()