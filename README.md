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
- Full support for DR and DQ loci serological nomenclature
- Broad/split antigen relationship handling
- 99.4% compatibility with laboratory bead mapping files (160/161 non-DP mappings)

## Installation

### Install from GitHub

```bash
# Install latest version from GitHub
pip install git+https://github.com/tavinathanson/imgtsero.git

# Install a specific version (recommended for reproducibility)
pip install git+https://github.com/tavinathanson/imgtsero.git@v0.3.1

# Or install from main branch
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

Download IMGT/HLA data files for a specific version (optional - data is auto-downloaded when needed):

```bash
# Using the installed command
imgtsero 3610

# Or using Python module
python -m imgtsero 3610
```

This downloads:
- `data/rel_dna_ser.3610.txt` - WMDA molecular to serological mapping for version 3.61.0
- `data/rel_ser_ser.3610.txt` - WMDA serological antigen relationships for version 3.61.0

**Note**: You no longer need to manually download data - it's automatically downloaded when you initialize the converter with a version.

### Python API

```python
import imgtsero

# Initialize with IMGT/HLA database version (data auto-downloaded if needed)
converter = imgtsero.HLAConverter(3610)  # Version 3.61.0

# Or use convenience function (version required)
result = imgtsero.convert("A*01:01", 3610)
# Returns: "A1" (molecular -> serological)

# Bidirectional HLA conversion
# Auto-detect format and convert to opposite
result = converter.convert("A*01:01")
# Returns: "A1" (molecular -> serological)

result = converter.convert("A1")
# Returns: ["A*01:01", "A*01:02", "A*01:03"] (serological -> molecular)

# Explicit format conversion
result = converter.convert("A*01:01", "s")  # Convert to serological
# Returns: "A1"

result = converter.convert("A1", "m")  # Convert to molecular
# Returns: ["A*01:01", "A*01:02", "A*01:03"]

# Broad/split antigen support
# Normal: A*02:03 maps to split antigen A203
result = converter.convert("A*02:03")
# Returns: "A203"

# With return_broad=True: return broad antigen A2 instead
result = converter.convert("A*02:03", return_broad=True)
# Returns: "A2"

# Normal: A2 maps to direct A2 alleles only
result = converter.convert("A2")
# Returns: ["A*02:01", "A*02:02", ...] (56 alleles)

# With expand_splits=True: include split antigens A203, A210
result = converter.convert("A2", expand_splits=True)
# Returns: ["A*02:01", "A*02:02", "A*02:03", "A*02:10", ...] (58 alleles)

# Using convenience function with version
result = imgtsero.convert("A2", 3610, expand_splits=True)
# Returns: ["A*02:01", "A*02:02", "A*02:03", "A*02:10", ...] (58 alleles)

# Error handling - raises HLAConversionError for unrecognized alleles
try:
    result = converter.convert("A*99:99")  # Invalid allele
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

## Development and Release Process

### Running Tests

```bash
# Run all tests
python3 -m unittest discover tests/

# Run specific test modules
python3 -m unittest tests.test_converter -v
python3 -m unittest tests.test_parser -v
```

### Creating Releases

To create a new release version:

1. **Update version numbers** in all relevant files:
   ```bash
   # Update version in imgtsero/__init__.py
   __version__ = "0.3.1"
   
   # Update version in pyproject.toml
   version = "0.3.1"
   
   # Update version in setup.py
   version="0.3.1"
   ```

2. **Commit the version changes**:
   ```bash
   git add imgtsero/__init__.py pyproject.toml setup.py
   git commit -m "Bump version to 0.3.1"
   ```

3. **Create and push a git tag**:
   ```bash
   # Create an annotated tag
   git tag -a v0.3.1 -m "Release v0.3.1: Add DR/DQ support and bead mapping compatibility"
   
   # Push the tag to origin
   git push origin v0.3.1
   
   # Or push all tags
   git push --tags
   ```

4. **Verify the release**:
   ```bash
   # Check that the tag was created
   git tag -l
   
   # Install from the specific tag to test
   pip install git+https://github.com/tavinathanson/imgtsero.git@v0.3.1
   ```

### Installation from Specific Versions

Users can install specific versions using git tags:

```bash
# Latest release
pip install git+https://github.com/tavinathanson/imgtsero.git@v0.3.1

# Previous versions
pip install git+https://github.com/tavinathanson/imgtsero.git@v0.3.0
pip install git+https://github.com/tavinathanson/imgtsero.git@v0.2.1
pip install git+https://github.com/tavinathanson/imgtsero.git@v0.2.0

# Development branch
pip install git+https://github.com/tavinathanson/imgtsero.git@main
```

## Implementation Notes

### HLA Serological Nomenclature Standards

This library follows **official WHO/IMGT nomenclature standards** as implemented in the WMDA data files. We do not invent any custom naming schemes.

#### C Locus Nomenclature

The C locus requires special handling due to historical naming conventions:

- **Serological format**: Uses "Cw" prefix (e.g., Cw1, Cw14, Cw12) to avoid confusion with complement components
- **Molecular format**: Uses "C*" prefix without the "w" (e.g., C*01:02, C*14:02, C*12:02)

**Implementation based on WMDA data:**
- C locus alleles with explicit serological numbers (e.g., "1", "2") → Cw1, Cw2
- C locus alleles with "?" in serological field → Use first allele field (C*14:02 → Cw14)

#### DR Locus Nomenclature

DR serological antigens follow **official WHO nomenclature**:

- **DRB1** alleles → **DR1-DR17** serological antigens (e.g., DRB1*01:01 → DR1)
- **DRB3** alleles → **DR52** serological antigen  
- **DRB4** alleles → **DR53** serological antigen
- **DRB5** alleles → **DR51** serological antigen

**Examples:**
```python
converter.convert("DRB1*01:01")  # Returns "DR1"
converter.convert("DRB1*04:01")  # Returns "DR4" 
converter.convert("DRB5*01:01")  # Returns "DR51"
converter.convert("DR1")         # Returns ["DRB1*01:01", "DRB1*01:02", ...]
```

#### DQ Locus Nomenclature

DQ serological antigens correspond to **DQB1** alleles in WMDA data:

- **DQB1** alleles → **DQ1-DQ9** serological antigens (e.g., DQB1*02:01 → DQ2)

**Note:** While actual HLA-DQ molecules are heterodimers requiring both DQA1 and DQB1 chains, the WMDA rel_dna_ser.txt file only provides serological equivalents for DQB1 alleles, following standard practice.

#### DP Locus Limitation

**DP antigens are not supported** because:
- WMDA data shows DPB1 alleles have "?" for serological equivalents
- DP molecules are heterodimers requiring both DPA1 and DPB1 chains
- No standardized serological nomenclature exists in the WMDA data

**Standards Compliance:**
- All nomenclature follows WHO Nomenclature Committee for Factors of the HLA System
- Data sourced from official WMDA files maintained by IMGT/HLA Database  
- No custom or invented naming schemes are used
- Parser correctly implements WMDA rel_dna_ser.txt column format:
  - Column 3: Unambiguous Serological Antigen (preferred)
  - Column 5: Assumed Serological Antigen (fallback)  
  - Column 6: Expert assigned exceptions (final fallback)
- Priority order: Unambiguous > Assumed > Expert (per WMDA documentation)

**Bead Mapping Compatibility:**
- 100% compatibility with WMDA-standard bead mappings  
- Some laboratory bead mappings may differ from WMDA standards (e.g., assay-specific allele selections)
- Heterodimer specifications (DQA1+DQB1) not supported as WMDA provides only single-chain mappings

### Broad/Split Antigen Support

The library supports broad and split antigen relationships as defined by WMDA:

**Broad antigens** are serological specificities that encompass multiple split antigens:
- Example: A2 is a broad antigen that includes splits A203 and A210

**Split antigens** are more specific serological types within a broad antigen:
- Example: A203 and A210 are splits of the broad antigen A2

**Parameters:**
- `expand_splits=True`: When converting from broad antigen to molecular, include alleles from all split antigens
- `return_broad=True`: When converting from molecular to serological, return the broad antigen instead of split

**Examples:**
```python
# Without broad/split support
imgtsero.convert("A2")           # Returns only direct A2 alleles
imgtsero.convert("A*02:03")      # Returns "A203" (split)

# With broad/split support  
imgtsero.convert("A2", expand_splits=True)      # Includes A203 and A210 alleles
imgtsero.convert("A*02:03", return_broad=True)  # Returns "A2" (broad)
```

### Data Sources

The library uses official WMDA (World Marrow Donor Association) files:
- `rel_dna_ser.txt`: Maps molecular alleles to serological equivalents
- `rel_ser_ser.txt`: Maps broad and split serological relationships

### Bead Mapping Compatibility

The library has been extensively tested against laboratory bead mapping files and achieves:
- **99.4% compatibility** (160/161) for all non-DP mappings
- Full support for DR locus mappings (DR1-DR17, DR51-DR53)
- Full support for DQ locus mappings (DQ1-DQ9)
- Automatic extraction of DQB1 alleles from heterodimer specifications
- The single incompatible mapping (DQ7 → DQB1*03:19) is due to an error in the bead mapping file where DQB1*03:19 actually maps to DQ3 in WMDA data

## License

MIT License