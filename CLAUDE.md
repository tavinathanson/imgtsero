# Claude Code Instructions for imgtsero

## Data Storage

All cached data (both WMDA files and KIR ligand cache) is stored in the `data_dir` parameter (default: "data"):
- WMDA files: `rel_dna_ser.{version}.txt` and `rel_ser_ser.{version}.txt`
- KIR ligand cache: `kir_ligand_{version}.json`

## Testing

This project has a flake.nix file. Always use `nix develop` for testing:

```bash
# First, enter the development shell
nix develop

# Then run tests normally
pytest tests/ -v
```

Examples in nix develop shell:
```bash
# Run all tests
pytest tests/ -v

# Run specific test file
pytest tests/test_kir_ligand.py -v

# Run specific test class
pytest tests/test_kir_ligand.py::TestKIRLigandCompression -v

# Run specific test
pytest tests/test_kir_ligand.py::TestKIRLigandCompression::test_compression_adds_four_digit_forms -v

# Run the new compression tests
pytest tests/test_kir_ligand.py::TestKIRLigandCompression::test_compression_skips_none_values -v
pytest tests/test_kir_ligand.py::TestKIRLigandCompression::test_c17_01_classification_with_real_data -v
pytest tests/test_kir_ligand.py::TestKIRLigandCompression::test_api_data_completeness -v
```

## Development Environment

- This project uses `nix develop` which automatically creates and activates a .venv
- The flake.nix ensures consistent Python environment and installs requirements.txt
- Do NOT use plain `nix-shell` - use `nix develop` instead

## Recent Updates (v0.4.8)

### KIR Ligand Classification Improvements

1. **Pagination Support**: The API fetch now handles pagination to retrieve all ~28,000+ alleles instead of just the first 1,000
2. **Null Value Handling**: The compression logic now skips None/null values when checking consistency
   - Example: C*17:01 now correctly compresses to C2 even though C*17:01:01:01 has null KIR data
3. **Data Completeness**: Added tests to ensure we get the expected number of alleles from the API

### Testing KIR Classification

To test a specific allele's KIR classification:
```python
import imgtsero
converter = imgtsero.HLAConverter(version=3610, enable_kir=True)
result = converter.classify_kir_ligand("C*17:01")
print(f"KIR type: {result['kir_ligand_type']}")  # Should print "C2"
```