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