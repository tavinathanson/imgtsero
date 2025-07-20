"""HLA format conversion functionality."""

import re
from typing import Union, List, Optional
from .parser import HLAParser


class HLAConversionError(Exception):
    """Raised when HLA conversion fails due to unrecognized allele."""
    pass


class HLAConverter:
    """Bidirectional HLA format converter."""
    
    def __init__(self, version: int, data_dir="data"):
        self.parser = HLAParser(version, data_dir)
    
    def convert(self, hla_type: str, target_format: Optional[str] = None, 
                expand_splits: bool = False, return_broad: bool = False) -> Union[str, List[str]]:
        """
        Convert HLA typing between serological and molecular formats.
        
        Args:
            hla_type: Input HLA type (e.g., "A*01:01", "A1")
            target_format: Target format ("s" for serological, "m" for molecular, 
                          or None for auto-detect and convert to opposite)
            expand_splits: When converting serological to molecular, include alleles 
                          from split antigens if this is a broad antigen (e.g., A2 includes A203)
            return_broad: When converting molecular to serological, return broad 
                         antigen instead of split if available (e.g., A*02:03 -> A2 instead of A203)
        
        Returns:
            Converted HLA type(s). Returns string for serological, list for molecular.
        """
        if target_format is None:
            # Auto-detect and convert to opposite
            if self._is_serological(hla_type):
                return self._to_molecular(hla_type, expand_splits)
            else:
                return self._to_serology(hla_type, return_broad)
        elif target_format == "s":
            return self._to_serology(hla_type, return_broad)
        elif target_format == "m":
            return self._to_molecular(hla_type, expand_splits)
        else:
            raise ValueError(f"Unsupported target format: {target_format}. Use 's' for serological, 'm' for molecular, or None for auto-detect.")
    
    def _to_2field(self, molecular_allele: str) -> str:
        """Convert molecular allele to 2-field format (e.g., A*01:01:01:01 -> A*01:01)."""
        if '*' in molecular_allele and ':' in molecular_allele:
            parts = molecular_allele.split(':')
            if len(parts) >= 2:
                return f"{parts[0]}:{parts[1]}"
        
        # If it's not a molecular allele or doesn't have enough fields, return as-is
        return molecular_allele
    
    def _to_serology(self, hla_type: str, return_broad: bool = False) -> str:
        """Convert molecular allele to serological equivalent."""
        if '*' in hla_type:
            # Validate that this is a real molecular allele first
            if not self._is_valid_molecular_allele(hla_type):
                raise HLAConversionError(f"Unrecognized molecular allele: {hla_type}")
            
            # Try exact match first
            sero = self.parser.get_molecular_to_serological(hla_type, return_broad)
            if sero:
                return sero
            
            # Try 2-field version
            two_field = self._to_2field(hla_type)
            sero = self.parser.get_molecular_to_serological(two_field, return_broad)
            if sero:
                return sero
            
            # If we have a valid molecular allele but no serological mapping
            raise HLAConversionError(f"No serological equivalent found for molecular allele: {hla_type}")
        
        # If already serological, validate and return as-is
        if self._is_serological(hla_type):
            if not self._is_valid_serological_allele(hla_type):
                raise HLAConversionError(f"Unrecognized serological allele: {hla_type}")
            return hla_type
        
        # Invalid format
        raise HLAConversionError(f"Invalid HLA format: {hla_type}")    
    
    def _to_molecular(self, hla_type: str, expand_splits: bool = False) -> List[str]:
        """Convert serological to molecular alleles (2-field format)."""
        if self._is_serological(hla_type):
            # Validate serological allele
            if not self._is_valid_serological_allele(hla_type):
                raise HLAConversionError(f"Unrecognized serological allele: {hla_type}")
            
            # Convert all to 2-field format
            raw_alleles = self.parser.get_serological_mapping(hla_type, expand_splits)
            if not raw_alleles:
                raise HLAConversionError(f"No molecular equivalents found for serological allele: {hla_type}")
            return [self._to_2field(allele) for allele in raw_alleles]
        
        # If already molecular, validate and return as 2-field in list
        if '*' in hla_type:
            if not self._is_valid_molecular_allele(hla_type):
                raise HLAConversionError(f"Unrecognized molecular allele: {hla_type}")
            return [self._to_2field(hla_type)]
        
        # Invalid format
        raise HLAConversionError(f"Invalid HLA format: {hla_type}")
    
    def _is_serological(self, hla_type: str) -> bool:
        """Check if HLA type is in serological format."""
        # Serological types typically don't contain '*' and may contain numbers
        # Handle standard format (A1, B27) and C locus format (Cw14)
        return '*' not in hla_type and (re.match(r'^[A-Z]+\d+$', hla_type) or re.match(r'^Cw\d+$', hla_type))
    
    def _is_valid_molecular_allele(self, hla_type: str) -> bool:
        """Check if molecular allele exists in the database."""
        if not '*' in hla_type:
            return False
        
        # Extract locus
        locus_match = re.match(r'^([A-Z]+\d*)\*', hla_type)
        if not locus_match:
            return False
        locus = locus_match.group(1)
        
        # Check if allele exists in the database
        alleles = self.parser.get_alleles_for_locus(locus)
        
        # Check exact match
        if hla_type in alleles:
            return True
        
        # Check if 2-field version exists
        two_field = self._to_2field(hla_type)
        for allele in alleles:
            if self._to_2field(allele) == two_field:
                return True
        
        return False
    
    def _is_valid_serological_allele(self, hla_type: str) -> bool:
        """Check if serological allele is known."""
        # Ensure data is loaded
        self.parser._load_data()
        # Check if it's in our serological mapping
        return hla_type in self.parser._serological_mapping


def convert(hla_type: str, version: int, target_format: Optional[str] = None, data_dir: str = "data",
            expand_splits: bool = False, return_broad: bool = False) -> Union[str, List[str]]:
    """
    Convenience function for HLA conversion.
    
    Args:
        hla_type: Input HLA type (e.g., "A*01:01", "A1")
        version: IMGT/HLA database version (e.g., 3610 for version 3.61.0)
        target_format: Target format ("s" for serological, "m" for molecular, 
                      or None for auto-detect and convert to opposite)
        data_dir: Directory containing HLA data files
        expand_splits: When converting serological to molecular, include alleles 
                      from split antigens if this is a broad antigen (e.g., A2 includes A203)
        return_broad: When converting molecular to serological, return broad 
                     antigen instead of split if available (e.g., A*02:03 -> A2 instead of A203)
    
    Returns:
        Converted HLA type(s). Returns string for serological, list for molecular.
    """
    converter = HLAConverter(version, data_dir)
    return converter.convert(hla_type, target_format, expand_splits, return_broad)