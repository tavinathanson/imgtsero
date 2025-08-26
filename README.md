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
- KIR ligand classification for HLA alleles using IPD-IMGT/HLA API data

## Installation

### Install from GitHub

```bash
# Install latest version from GitHub
pip install git+https://github.com/tavinathanson/imgtsero.git

# Install a specific version (recommended for reproducibility)
pip install git+https://github.com/tavinathanson/imgtsero.git@v0.4.0

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
# Default: A*02:03 maps to split antigen A203
result = converter.convert("A*02:03")
# Returns: "A203"

# Display options for broad/split relationships
result = converter.convert("A*02:03", handle_broad="broad")
# Returns: "A2" (show broad only)

result = converter.convert("A*02:03", handle_broad="both")  
# Returns: "A2 (A203)" (show both broad and split)

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

### KIR Ligand Classification

The library can classify HLA alleles as KIR (Killer-cell Immunoglobulin-like Receptor) ligands using data from the IPD-IMGT/HLA API.

#### KIR Ligand Biology

- **HLA-B**: All alleles have either Bw4 or Bw6 epitope
  - Bw4 is a KIR ligand (binds KIR3DL1)
  - Bw6 is NOT a KIR ligand
- **HLA-C**: All alleles are either C1 or C2 group
  - C1 binds KIR2DL2/2DL3
  - C2 binds KIR2DL1
  - Both are KIR ligands
- **HLA-A**: Most are not KIR ligands, except some with Bw4 (A*23, A*24, A*32)

#### Basic Usage

```python
from imgtsero import HLAConverter

# Initialize with KIR classification enabled
# Version can be either format: 3610 or "3.61.0"
converter = HLAConverter(version=3610, enable_kir=True)

# Classify a molecular allele
result = converter.classify_kir_ligand("B*27:05:02")
print(result)
# {
#     "is_kir_ligand": True,
#     "kir_ligand_type": "Bw4",
#     "kir_receptors": ["KIR3DL1"],
#     "source": "api"
# }

# Classify a serological antigen
result = converter.classify_kir_ligand("Cw7")
# Returns: {"is_kir_ligand": True, "kir_ligand_type": "C2", ...}

# Parse and validate SAB bead annotations
from imgtsero import KIRLigandClassifier

bead_name = "B27,Bw4"
hla_antigen, kir_annotation = KIRLigandClassifier.parse_bead_annotation(bead_name)
# Returns: ("B27", "Bw4")

# Classify with validation - raises ValueError if annotation conflicts with API
try:
    result = converter.classify_kir_ligand(hla_antigen, bead_annotation=kir_annotation)
except ValueError as e:
    print(f"Data conflict: {e}")
```

#### Direct KIR Classifier Usage

```python
from imgtsero import KIRLigandClassifier

# Initialize classifier (handles version format conversion)
classifier = KIRLigandClassifier(version=3610)  # or "3.61.0"
classifier.load_data()

# Get all alleles grouped by KIR ligand type
kir_groups = classifier.get_all_kir_ligands()
# Returns: {
#     "Bw4": ["B*27:05", "B*44:02", "A*23:01", ...],
#     "Bw6": ["B*07:02", "B*08:01", ...],
#     "C1": ["C*01:02", "C*03:03", ...],
#     "C2": ["C*02:02", "C*04:01", ...]
# }
```

#### Implementation Details

- **Data Source**: IPD-IMGT/HLA API provides official KIR ligand assignments
- **Automatic pagination**: The API fetch handles pagination to retrieve all HLA-A, B, and C alleles (~28,222 alleles for version 3.61.0) across multiple requests. Note: This represents only A*, B*, and C* loci out of the total 43,225+ alleles in the full database (which includes DRB1*, DQB1*, and many other loci)
- **Automatic 4-digit compression**: When the API returns only high-resolution alleles (e.g., B*27:05:02, B*27:05:09), the system automatically creates 4-digit entries (e.g., B*27:05) if all high-resolution variants with non-null KIR data have consistent KIR ligand types
- **Null value handling**: Alleles with no KIR ligand data (null/None values) are ignored during consistency checking. For example, if C*17:01:01:01 has no data but C*17:01:01:02 and C*17:01:02 both have "C2", then C*17:01 will be compressed to "C2"
- **Consistency validation**: If high-resolution alleles of the same 4-digit type have conflicting KIR ligand assignments (excluding null values), a ValueError is raised with detailed information about the conflicts
- **Caching**: Data is cached locally after first retrieval for performance (includes all paginated results)
- **Validation**: SAB bead annotations (e.g., "Bw4", "Bw6") are used only for validation against API data
- **Error Handling**: Raises ValueError if bead annotation conflicts with API data or if inconsistent KIR ligand types are found during compression
- **Version Support**: Handles both 3610 and "3.61.0" version formats

## Testing

```bash
# Run all tests (excluding slow tests)
pytest

# Run all tests including slow tests  
pytest -m ""

# Run only the slow pagination test (fetches all 43,416 HLA alleles)
pytest -m slow

# Run a specific test
pytest tests/test_kir_ligand.py::TestKIRLigandCompression::test_c17_01_classification_with_real_data
```

The slow test `test_total_database_alleles_by_pagination_v3_61_0` fetches all 43,416 HLA alleles across 44 API pages to verify our pagination implementation. It takes about 40 seconds to complete.

## Project Structure

```
imgtsero/
├── imgtsero/
│   ├── __init__.py          # Main API exports
│   ├── __main__.py          # CLI entry point
│   ├── downloader.py        # Data download functionality
│   ├── parser.py            # HLA data parsing
│   ├── converter.py         # HLA format conversion
│   └── kir_ligand.py        # KIR ligand classification
├── tests/
│   ├── test_downloader.py
│   ├── test_parser.py
│   ├── test_converter.py
│   └── test_kir_ligand.py
├── README.md
├── setup.py
└── requirements.txt
```

## Requirements

- Python 3.7+
- No external dependencies (uses only Python standard library)

## Testing

### Using nix develop (recommended if you have flake.nix)

```bash
# Enter the development shell (will create/activate .venv automatically)
nix develop

# Then run tests
pytest tests/ -v

# Or run specific tests
pytest tests/test_kir_ligand.py -v
pytest tests/test_kir_ligand.py::TestKIRLigandCompression -v
```

### Using nix-shell (without flake)

```bash
# Run all tests
nix-shell -p python3 python3Packages.pytest --run "python -m pytest tests/ -v"

# Run specific test file
nix-shell -p python3 python3Packages.pytest --run "python -m pytest tests/test_kir_ligand.py -v"
```

### Using local environment

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

1. **Run the version bump script**:
   ```bash
   # Bump version (updates all files, commits, and tags)
   ./bump_version.py 0.4.1
   
   # Or with a custom tag message
   ./bump_version.py 0.4.1 -m "Release v{version}: Add new feature"
   
   # Or update files only (no git operations)
   ./bump_version.py 0.4.1 --no-git
   ```

2. **Push to GitHub**:
   ```bash
   git push origin main
   git push origin v0.4.1
   ```

3. **Verify the release**:
   ```bash
   # Install from the specific tag to test
   pip install git+https://github.com/tavinathanson/imgtsero.git@v0.4.1
   ```

### Installation from Specific Versions

Users can install specific versions using git tags:

```bash
# Latest release
pip install git+https://github.com/tavinathanson/imgtsero.git@v0.4.0

# Previous versions
pip install git+https://github.com/tavinathanson/imgtsero.git@v0.3.2
pip install git+https://github.com/tavinathanson/imgtsero.git@v0.3.1
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
- `handle_broad`: Controls molecular to serological display format:
  - `"split"` (default): Show split antigen (e.g., "A203")
  - `"broad"`: Show broad antigen when available (e.g., "A2")
  - `"both"`: Show "Broad (Split)" format when both exist (e.g., "A2 (A203)")

**Examples:**
```python
# Serological to molecular conversion
imgtsero.convert("A2")                          # Returns only direct A2 alleles
imgtsero.convert("A2", expand_splits=True)      # Includes A203 and A210 alleles

# Molecular to serological conversion options
imgtsero.convert("A*02:03")                     # Returns "A203" (split, default)
imgtsero.convert("A*02:03", handle_broad="broad")  # Returns "A2" (broad only)
imgtsero.convert("A*02:03", handle_broad="both")   # Returns "A2 (A203)" (combined)
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
- The single incompatible mapping is a confirmed error in the bead mapping file:
  - **Bead mapping claims**: DQ7 → DQB1*03:19 (line 155: `DQ7,"DQA1*05:05, DQB1*03:19"`)
  - **WMDA data shows**: DQB1*03:19 maps to DQ3 (has serological "3" in column 5)
  - **Correct DQ7 mappings**: DQB1*03:01 and DQB1*03:04 (have serological "7" in column 3)
  - This is a data error where the bead file incorrectly associates DQB1*03:19 with DQ7

## License

MIT License