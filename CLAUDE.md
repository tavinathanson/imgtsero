# Claude Code Instructions for imgtsero

## Testing

This project uses nix-shell for testing. Always run tests using:

```bash
nix-shell -p python3 python3Packages.pytest --run "python -m pytest [test_file_or_directory] -v"
```

Examples:
```bash
# Run all tests
nix-shell -p python3 python3Packages.pytest --run "python -m pytest tests/ -v"

# Run specific test file
nix-shell -p python3 python3Packages.pytest --run "python -m pytest tests/test_kir_ligand.py -v"

# Run specific test class
nix-shell -p python3 python3Packages.pytest --run "python -m pytest tests/test_kir_ligand.py::TestKIRLigandCompression -v"

# Run specific test
nix-shell -p python3 python3Packages.pytest --run "python -m pytest tests/test_kir_ligand.py::TestKIRLigandCompression::test_compression_adds_four_digit_forms -v"
```

## Development Environment

Do NOT use plain `python` or `pytest` commands - always use nix-shell as shown above.