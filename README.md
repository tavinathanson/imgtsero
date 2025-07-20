# imgtsero

A Python library for downloading and working with IMGT/HLA data files from the official GitHub repository.

**Warning**: not yet intended for general use. May have bugs: please inquire if interested!

## Features

- Download IMGT/HLA WMDA files from official repository for comprehensive serological mapping
- Bidirectional conversion between serological and molecular HLA formats
- Auto-detection of input format with conversion to opposite
- Robust error handling for unrecognized alleles
- Command-line interface for data downloads
- Uses official WHO/WMDA serological mapping data for accurate conversions
- Supports 26+ HLA loci with thousands of allele-to-serotype mappings
- Special handling for C locus nomenclature (Cw format)

## Installation

### Install from GitHub

```bash
# Install directly from GitHub
pip install git+https://github.com/tavinathanson/imgtsero.git

# Or install a specific branch/tag
pip install git+https://github.com/tavinathanson/imgtsero.git@main
```

### Development Installation

```bash
# Clone the repository
git clone https://github.com/tavinathanson/imgtsero.git
cd imgtsero

# Install in development mode with dev dependencies
pip install -e ".[dev]"

# Or just install dev dependencies
pip install -r requirements.txt
```

## Usage

### Command Line Interface

After installation, download IMGT/HLA data files for a specific version:

```bash
# Using the installed command
imgtsero 3610

# Or using Python module
python -m imgtsero 3610
```

This downloads:
- `data/rel_dna_ser.3610.txt` - WMDA molecular to serological mapping for version 3.61.0
- `data/rel_ser_ser.3610.txt` - WMDA serological antigen relationships for version 3.61.0

### Python API

```python
import imgtsero

# Download data files
imgtsero.download_data(3610)

# Bidirectional HLA conversion
# Auto-detect format and convert to opposite
result = imgtsero.convert("A*01:01")
# Returns: "A1" (molecular -> serological)

result = imgtsero.convert("A1")
# Returns: ["A*01:01", "A*01:02", "A*01:03"] (serological -> molecular)

# Explicit format conversion
result = imgtsero.convert("A*01:01", "s")  # Convert to serological
# Returns: "A1"

result = imgtsero.convert("A1", "m")  # Convert to molecular
# Returns: ["A*01:01", "A*01:02", "A*01:03"]

# Broad/split antigen support
# Normal: A*02:03 maps to split antigen A203
result = imgtsero.convert("A*02:03")
# Returns: "A203"

# With prefer_broad=True: return broad antigen A2 instead
result = imgtsero.convert("A*02:03", prefer_broad=True)
# Returns: "A2"

# Normal: A2 maps to direct A2 alleles only
result = imgtsero.convert("A2")
# Returns: ["A*02:01", "A*02:02", ...] (56 alleles)

# With include_broad=True: include split antigens A203, A210
result = imgtsero.convert("A2", include_broad=True)
# Returns: ["A*02:01", "A*02:02", "A*02:03", "A*02:10", ...] (58 alleles)

# Error handling - raises HLAConversionError for unrecognized alleles
try:
    result = imgtsero.convert("A*99:99")  # Invalid allele
except imgtsero.HLAConversionError as e:
    print(f"Error: {e}")  # "Unrecognized molecular allele: A*99:99"
```

## Project Structure

```
imgtsero/
├── imgtsero/
│   ├── __init__.py          # Main API exports
│   ├── __main__.py          # CLI entry point
│   ├── downloader.py        # Data download functionality
│   ├── parser.py            # HLA data parsing
│   └── converter.py         # HLA format conversion
├── tests/
│   ├── test_downloader.py
│   ├── test_parser.py
│   └── test_converter.py
├── README.md
├── setup.py
└── requirements.txt
```

## Requirements

- Python 3.7+
- No external dependencies (uses only Python standard library)

## Testing

```bash
python -m pytest tests/
```

## Implementation Notes

### C Locus Nomenclature

The C locus requires special handling due to historical naming conventions:

- **Serological format**: Uses "Cw" prefix (e.g., Cw14, Cw12) to avoid confusion with complement components
- **Molecular format**: Uses "C*" prefix without the "w" (e.g., C*14:02, C*12:02)

In the WMDA `rel_dna_ser.txt` file, C locus alleles have "?" in the serological column rather than numeric values. Our implementation handles this by:

1. Detecting C locus alleles with "?" in the serological field
2. Extracting the first field of the allele (e.g., "14" from C*14:02:01:01)
3. Mapping to the corresponding Cw serological name (e.g., Cw14)

This approach is consistent with:
- The WMDA `rel_ser_ser.txt` file which uses "Cw" nomenclature
- HLA nomenclature standards where the first field typically corresponds to the serological specificity
- Historical HLA workshop designations where "w" was retained for C locus antigens

### Broad/Split Antigen Support

The library supports broad and split antigen relationships as defined by WMDA:

**Broad antigens** are serological specificities that encompass multiple split antigens:
- Example: A2 is a broad antigen that includes splits A203 and A210

**Split antigens** are more specific serological types within a broad antigen:
- Example: A203 and A210 are splits of the broad antigen A2

**Parameters:**
- `include_broad=True`: When converting from broad antigen to molecular, include alleles from all split antigens
- `prefer_broad=True`: When converting from molecular to serological, return the broad antigen instead of split

**Examples:**
```python
# Without broad/split support
imgtsero.convert("A2")           # Returns only direct A2 alleles
imgtsero.convert("A*02:03")      # Returns "A203" (split)

# With broad/split support  
imgtsero.convert("A2", include_broad=True)      # Includes A203 and A210 alleles
imgtsero.convert("A*02:03", prefer_broad=True)  # Returns "A2" (broad)
```

### Data Sources

The library uses official WMDA (World Marrow Donor Association) files:
- `rel_dna_ser.txt`: Maps molecular alleles to serological equivalents
- `rel_ser_ser.txt`: Maps broad and split serological relationships

## License

MIT License