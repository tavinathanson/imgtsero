"""HLA format conversion functionality."""

import re
from typing import Union, List, Optional
from .parser import HLAParser


class HLAConversionError(Exception):
    """Raised when HLA conversion fails due to unrecognized allele."""
    pass


class HLAConverter:
    """Bidirectional HLA format converter."""
    
    def __init__(self, data_dir="data"):
        self.parser = HLAParser(data_dir)
    
    def convert(self, hla_type: str, target_format: Optional[str] = None) -> Union[str, List[str]]:
        """
        Convert HLA typing between serological and molecular formats.
        
        Args:
            hla_type: Input HLA type (e.g., "A*01:01", "A1")
            target_format: Target format ("s" for serological, "m" for molecular, 
                          or None for auto-detect and convert to opposite)
        
        Returns:
            Converted HLA type(s). Returns string for serological, list for molecular.
        """
        if target_format is None:
            # Auto-detect and convert to opposite
            if self._is_serological(hla_type):
                return self._to_molecular(hla_type)
            else:
                return self._to_serology(hla_type)
        elif target_format == "s":
            return self._to_serology(hla_type)
        elif target_format == "m":
            return self._to_molecular(hla_type)
        else:
            raise ValueError(f"Unsupported target format: {target_format}. Use 's' for serological, 'm' for molecular, or None for auto-detect.")
    
    def _to_2field(self, hla_type: str) -> str:
        """Convert to 2-field molecular format (e.g., A*01:01)."""
        # If already molecular, truncate to 2 fields
        if '*' in hla_type and ':' in hla_type:
            parts = hla_type.split(':')
            if len(parts) >= 2:
                return f"{parts[0]}:{parts[1]}"
        
        # If serological, convert to molecular first then to 2-field
        if self._is_serological(hla_type):
            molecular_alleles = self.parser.get_serological_mapping(hla_type)
            if molecular_alleles:
                # Return the first common allele in 2-field format
                return self._to_2field(molecular_alleles[0])
        
        return hla_type
    
    def _to_serology(self, hla_type: str) -> str:
        """Convert molecular allele to serological equivalent."""
        if '*' in hla_type:
            # Validate that this is a real molecular allele first
            if not self._is_valid_molecular_allele(hla_type):
                raise HLAConversionError(f"Unrecognized molecular allele: {hla_type}")
            
            # Try exact match first
            sero = self.parser.get_molecular_to_serological(hla_type)
            if sero:
                return sero
            
            # Try 2-field version
            two_field = self._to_2field(hla_type)
            sero = self.parser.get_molecular_to_serological(two_field)
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
    
    def _to_molecular(self, hla_type: str) -> List[str]:
        """Convert serological to molecular alleles (2-field format)."""
        if self._is_serological(hla_type):
            # Validate serological allele
            if not self._is_valid_serological_allele(hla_type):
                raise HLAConversionError(f"Unrecognized serological allele: {hla_type}")
            
            # Convert all to 2-field format
            raw_alleles = self.parser.get_serological_mapping(hla_type)
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


def convert(hla_type: str, target_format: Optional[str] = None, data_dir: str = "data") -> Union[str, List[str]]:
    """
    Convenience function for HLA conversion.
    
    Args:
        hla_type: Input HLA type (e.g., "A*01:01", "A1")
        target_format: Target format ("s" for serological, "m" for molecular, 
                      or None for auto-detect and convert to opposite)
        data_dir: Directory containing HLA data files
    
    Returns:
        Converted HLA type(s). Returns string for serological, list for molecular.
    """
    converter = HLAConverter(data_dir)
    return converter.convert(hla_type, target_format)