import imgtsero

def test_kir_ligand_detail():
    """Test that detailed KIR ligand subtypes are preserved"""
    converter = imgtsero.HLAConverter(version=3610, enable_kir=True)
    
    # Test cases showing Bw4 subtypes
    test_cases = [
        # Bw4 with 80T subtype
        ('B*27:05', 'Bw4', 'Bw4 - 80T'),
        ('B*44:02', 'Bw4', 'Bw4 - 80T'),
        ('B*44:03', 'Bw4', 'Bw4 - 80T'),
        
        # Bw4 with 80I subtype
        ('B*51:01', 'Bw4', 'Bw4 - 80I'),
        ('B*52:01', 'Bw4', 'Bw4 - 80I'),
        ('B*57:01', 'Bw4', 'Bw4 - 80I'),
        ('B*58:01', 'Bw4', 'Bw4 - 80I'),
        
        # Bw6 (no subtype)
        ('B*07:02', 'Bw6', 'Bw6'),
        
        # C ligands (no subtype)
        ('C*01:02', 'C1', 'C1'),
        ('C*04:01', 'C2', 'C2'),
    ]
    
    print("Testing KIR ligand classification with detail preservation:\n")
    for allele, expected_type, expected_detail in test_cases:
        result = converter.classify_kir_ligand(allele)
        print(f"{allele}:")
        print(f"  Base type: {result['kir_ligand_type']} (expected: {expected_type})")
        print(f"  Detailed type: {result['kir_ligand_type_detail']} (expected: {expected_detail})")
        print(f"  Is KIR ligand: {result['is_kir_ligand']}")
        
        # Verify results
        assert result['kir_ligand_type'] == expected_type, f"Base type mismatch for {allele}"
        assert result['kir_ligand_type_detail'] == expected_detail, f"Detail type mismatch for {allele}"
        
        # Verify is_kir_ligand is correct
        if expected_type in ['Bw4', 'C1', 'C2']:
            assert result['is_kir_ligand'] == True
        else:
            assert result['is_kir_ligand'] == False
        
        print()
    
    print("All detail preservation tests passed!")

if __name__ == '__main__':
    test_kir_ligand_detail()