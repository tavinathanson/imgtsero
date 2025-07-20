"""Parser for IMGT/HLA WMDA data files."""

import os
import re
from typing import Dict, List, Optional, Set


class HLAParser:
    """Parser for IMGT/HLA WMDA data files."""
    
    def __init__(self, data_dir="data"):
        self.data_dir = data_dir
        self._molecular_to_serological = {}  # Maps A*01:01 -> 1
        self._serological_to_molecular = {}  # Maps A1 -> [A*01:01, A*01:02, ...]
        self._serological_mapping = {}  # For backward compatibility
        self._allele_data = {}  # Maps locus -> set of alleles
        self._broad_to_splits = {}  # Maps A2 -> [A203, A210]
        self._split_to_broad = {}  # Maps A203 -> A2
        self._loaded = False
    
    def _load_data(self):
        """Load HLA data from WMDA files."""
        if self._loaded:
            return
        
        # Find the most recent WMDA files
        rel_dna_ser_files = [f for f in os.listdir(self.data_dir) if f.startswith("rel_dna_ser.") and f.endswith(".txt")]
        if not rel_dna_ser_files:
            raise FileNotFoundError("No rel_dna_ser files found in data directory")
        
        # Sort by version number (assume format rel_dna_ser.XXXX.txt)
        rel_dna_ser_files.sort(key=lambda x: int(re.search(r'rel_dna_ser\.(\d+)\.txt', x).group(1)), reverse=True)
        rel_dna_ser_file = os.path.join(self.data_dir, rel_dna_ser_files[0])
        
        # Load molecular to serological mapping
        self._load_rel_dna_ser(rel_dna_ser_file)
        
        # Load broad/split antigen relationships
        rel_ser_ser_files = [f for f in os.listdir(self.data_dir) if f.startswith("rel_ser_ser.") and f.endswith(".txt")]
        if rel_ser_ser_files:
            rel_ser_ser_files.sort(key=lambda x: int(re.search(r'rel_ser_ser\.(\d+)\.txt', x).group(1)), reverse=True)
            rel_ser_ser_file = os.path.join(self.data_dir, rel_ser_ser_files[0])
            self._load_rel_ser_ser(rel_ser_ser_file)
        
        self._loaded = True
    
    def _load_rel_dna_ser(self, filepath):
        """Load molecular to serological mapping from rel_dna_ser.txt."""
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    # Format: A*;01:01:01:01;1;;;
                    parts = line.split(';')
                    if len(parts) >= 3:
                        locus_prefix = parts[0]  # e.g., "A*", "C*"
                        allele_suffix = parts[1]  # e.g., "01:01:01:01"
                        serological_num = parts[2]  # e.g., "1" or "?" for C locus
                        
                        if locus_prefix and allele_suffix:
                            # Construct full allele name: A*01:01:01:01
                            full_allele = f"{locus_prefix}{allele_suffix}"
                            
                            # Extract locus for allele tracking
                            locus_match = re.match(r'^([A-Z]+\d*)\*', full_allele)
                            if locus_match:
                                locus = locus_match.group(1)
                                
                                # Track all alleles by locus
                                if locus not in self._allele_data:
                                    self._allele_data[locus] = set()
                                self._allele_data[locus].add(full_allele)
                                
                                # Convert to 2-field format for mapping
                                two_field = self._to_2field(full_allele)
                                
                                # Handle serological mapping
                                if serological_num and serological_num != '0':
                                    if locus == 'C' and serological_num == '?':
                                        # Special handling for C locus: use allele number for serological
                                        # C*14:02 -> Cw14, C*12:02 -> Cw12
                                        allele_parts = allele_suffix.split(':')
                                        if len(allele_parts) >= 2:
                                            sero_num = allele_parts[0]  # First field (e.g., "14", "12")
                                            self._molecular_to_serological[two_field] = sero_num
                                            sero_name = f"Cw{sero_num}"
                                            
                                            # Map serological to molecular (build lists)
                                            if sero_name not in self._serological_to_molecular:
                                                self._serological_to_molecular[sero_name] = set()
                                            self._serological_to_molecular[sero_name].add(two_field)
                                    elif serological_num != '?':
                                        # Standard mapping for other loci
                                        self._molecular_to_serological[two_field] = serological_num
                                        
                                        # Build serological name (e.g., "A1", "B27")
                                        sero_name = f"{locus}{serological_num}"
                                        
                                        # Map serological to molecular (build lists)
                                        if sero_name not in self._serological_to_molecular:
                                            self._serological_to_molecular[sero_name] = set()
                                        self._serological_to_molecular[sero_name].add(two_field)
        
        # Convert sets to lists and update backward compatibility mapping
        for sero_name, allele_set in self._serological_to_molecular.items():
            allele_list = sorted(list(allele_set))
            self._serological_to_molecular[sero_name] = allele_list
            self._serological_mapping[sero_name] = allele_list
    
    def _load_rel_ser_ser(self, filepath):
        """Load broad/split antigen relationships from rel_ser_ser.txt."""
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    # Format: A;2;;203/210
                    parts = line.split(';')
                    if len(parts) >= 4:
                        locus = parts[0]  # e.g., "A"
                        broad_num = parts[1]  # e.g., "2"
                        splits_field2 = parts[2]  # e.g., "23/24" or ""
                        splits_field3 = parts[3]  # e.g., "203/210" or ""
                        
                        if locus and broad_num:
                            broad_name = f"{locus}{broad_num}"  # e.g., "A2"
                            
                            # Handle C locus special case
                            if locus == "C":
                                broad_name = f"Cw{broad_num}"
                            
                            # Parse splits from both fields
                            all_splits = []
                            if splits_field2:
                                all_splits.extend(splits_field2.split('/'))
                            if splits_field3:
                                all_splits.extend(splits_field3.split('/'))
                            
                            # Build split names and mappings
                            split_names = []
                            for split_num in all_splits:
                                if split_num:  # Skip empty strings
                                    if locus == "C":
                                        split_name = f"Cw{split_num}"
                                    else:
                                        split_name = f"{locus}{split_num}"
                                    split_names.append(split_name)
                                    self._split_to_broad[split_name] = broad_name
                            
                            if split_names:
                                self._broad_to_splits[broad_name] = split_names
    
    def _to_2field(self, molecular_allele):
        """Convert molecular allele to 2-field format (e.g., A*01:01:01:01 -> A*01:01)."""
        if '*' in molecular_allele and ':' in molecular_allele:
            parts = molecular_allele.split(':')
            if len(parts) >= 2:
                return f"{parts[0]}:{parts[1]}"
        return molecular_allele
    
    def get_alleles_for_locus(self, locus: str) -> Set[str]:
        """Get all alleles for a specific locus."""
        self._load_data()
        return self._allele_data.get(locus, set())
    
    def get_loci(self) -> List[str]:
        """Get all available HLA loci."""
        self._load_data()
        return list(self._allele_data.keys())
    
    def find_alleles_by_pattern(self, pattern: str) -> List[str]:
        """Find alleles matching a pattern."""
        self._load_data()
        results = []
        for locus_alleles in self._allele_data.values():
            for allele in locus_alleles:
                if re.search(pattern, allele, re.IGNORECASE):
                    results.append(allele)
        return sorted(results)
    
    
    def get_molecular_to_serological(self, allele: str, prefer_broad: bool = False) -> Optional[str]:
        """Get serological equivalent for a molecular allele."""
        self._load_data()
        # Convert to 2-field format first
        two_field = self._to_2field(allele)
        serological_num = self._molecular_to_serological.get(two_field)
        if serological_num:
            # Extract locus
            locus_match = re.match(r'^([A-Z]+\d*)\*', allele)
            if locus_match:
                locus = locus_match.group(1)
                if locus == 'C':
                    # C locus uses Cw nomenclature
                    sero_name = f"Cw{serological_num}"
                else:
                    sero_name = f"{locus}{serological_num}"
                
                # If prefer_broad is True, try to return the broad antigen
                if prefer_broad and sero_name in self._split_to_broad:
                    return self._split_to_broad[sero_name]
                
                return sero_name
        return None
    
    def get_serological_mapping(self, serological: str, include_broad: bool = False) -> List[str]:
        """Get molecular alleles for a serological type."""
        self._load_data()
        alleles = list(self._serological_mapping.get(serological, []))
        
        # If include_broad is True and this is a broad antigen, include split antigens
        if include_broad and serological in self._broad_to_splits:
            for split_antigen in self._broad_to_splits[serological]:
                split_alleles = self._serological_mapping.get(split_antigen, [])
                alleles.extend(split_alleles)
            # Remove duplicates and sort
            alleles = sorted(list(set(alleles)))
        
        return alleles
    
    def get_broad_antigen(self, serological: str) -> Optional[str]:
        """Get broad antigen for a split antigen."""
        self._load_data()
        return self._split_to_broad.get(serological)
    
    def get_split_antigens(self, broad_antigen: str) -> List[str]:
        """Get split antigens for a broad antigen."""
        self._load_data()
        return self._broad_to_splits.get(broad_antigen, [])