# Claude Code Instructions for imgtsero

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
```

## Development Environment

- This project uses `nix develop` which automatically creates and activates a .venv
- The flake.nix ensures consistent Python environment and installs requirements.txt
- Do NOT use plain `nix-shell` - use `nix develop` instead