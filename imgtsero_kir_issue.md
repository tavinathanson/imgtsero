# KIR Ligand Classification Issue in imgtsero

## Instructions for Claude in imgtsero repo

Please investigate and fix the KIR ligand classification issue in imgtsero. The `classify_kir_ligand` method is not correctly identifying known KIR ligand alleles.

## Steps to investigate:

1. First, run this test script to reproduce the issue:

```python
import imgtsero

def test_kir_ligand_classification():
    """Test that common KIR ligand alleles are correctly classified"""
    converter = imgtsero.HLAConverter(version=3610, enable_kir=True)
    
    # Test cases with expected results
    test_cases = [
        # HLA-B alleles
        ('B*27:05', True, 'Bw4'),  # Known Bw4 allele
        ('B*44:02', True, 'Bw4'),  # Known Bw4 allele
        ('B*07:02', False, 'Bw6'), # Known Bw6 allele (not a KIR ligand)
        ('B*08:01', False, 'Bw6'), # Known Bw6 allele (not a KIR ligand)
        
        # HLA-C alleles (all are KIR ligands)
        ('C*01:02', True, 'C1'),   # Known C1 allele
        ('C*03:04', True, 'C1'),   # Known C1 allele
        ('C*04:01', True, 'C2'),   # Known C2 allele
        ('C*05:01', True, 'C2'),   # Known C2 allele
    ]
    
    failures = []
    for allele, expected_is_ligand, expected_type in test_cases:
        result = converter.classify_kir_ligand(allele)
        print(f"\nAllele: {allele}")
        print(f"Result: {result}")
        print(f"Expected: is_kir_ligand={expected_is_ligand}, kir_ligand_type='{expected_type}'")
        
        # Check results
        if result['is_kir_ligand'] != expected_is_ligand:
            failures.append(f"{allele}: is_kir_ligand should be {expected_is_ligand}")
        
        if expected_is_ligand and result['kir_ligand_type'] != expected_type:
            failures.append(f"{allele}: kir_ligand_type should be '{expected_type}'")
    
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
    else:
        print("\nAll tests passed!")

if __name__ == '__main__':
    test_kir_ligand_classification()
```

2. Investigate the code path for `classify_kir_ligand` method
   - Check how it determines if an allele is a KIR ligand
   - Look for any data loading or API calls that might be failing
   - Check if the KIR ligand data is being properly initialized when `enable_kir=True`

3. Expected behavior:
   - **Bw4** alleles (like B*27:05, B*44:02) should return:
     - `is_kir_ligand: True`
     - `kir_ligand_type: 'Bw4'`
     - `kir_receptors: ['KIR3DL1']` or similar
   
   - **Bw6** alleles (like B*07:02, B*08:01) should return:
     - `is_kir_ligand: False`
     - `kir_ligand_type: 'Bw6'`
     - `kir_receptors: []`
   
   - **C1** alleles (like C*01:02, C*03:04) should return:
     - `is_kir_ligand: True`
     - `kir_ligand_type: 'C1'`
     - `kir_receptors: ['KIR2DL2', 'KIR2DL3']` or similar
   
   - **C2** alleles (like C*04:01, C*05:01) should return:
     - `is_kir_ligand: True`
     - `kir_ligand_type: 'C2'`
     - `kir_receptors: ['KIR2DL1']` or similar

4. Background on KIR ligands:
   - **Bw4**: Public epitope on some HLA-B alleles that serves as a ligand for KIR3DL1
   - **Bw6**: Public epitope on other HLA-B alleles that is NOT a KIR ligand
   - **C1**: KIR ligand group for HLA-C alleles with asparagine at position 80, binds KIR2DL2/3
   - **C2**: KIR ligand group for HLA-C alleles with lysine at position 80, binds KIR2DL1
   - All HLA-C alleles are KIR ligands (either C1 or C2)

5. Current behavior (incorrect):
   All test alleles return:
   ```python
   {
       'is_kir_ligand': False,
       'kir_ligand_type': None,
       'kir_receptors': [],
       'source': None
   }
   ```

Please fix the KIR ligand classification to properly identify these alleles according to their known KIR ligand status.