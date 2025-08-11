"""HLA format conversion functionality."""

import re
from .parser import HLAParser
from .kir_ligand import KIRLigandClassifier


class HLAConversionError(Exception):
    """Raised when HLA conversion fails due to unrecognized allele."""
    pass


class HLAConverter:
    """Bidirectional HLA format converter."""
    
    def __init__(self, version, data_dir="data", enable_kir=False):
        self.parser = HLAParser(version, data_dir)
        self.kir_classifier = None
        if enable_kir:
            # KIRLigandClassifier handles version format conversion internally
            self.kir_classifier = KIRLigandClassifier(version, data_dir)
    
    def convert(self, hla_type, target_format=None, expand_splits=False, handle_broad="split"):
        """
        Convert HLA typing between serological and molecular formats.
        
        Args:
            hla_type: Input HLA type (e.g., "A*01:01", "A1")
            target_format: Target format ("s" for serological, "m" for molecular, 
                          or None for auto-detect and convert to opposite)
            expand_splits: When converting serological to molecular, include alleles 
                          from split antigens if this is a broad antigen (e.g., A2 includes A203)
            handle_broad: How to handle broad/split relationships in molecular to serological conversion:
                         "split" (default): Return split antigen (e.g., A*02:03 -> "A203")
                         "broad": Return broad antigen when available (e.g., A*02:03 -> "A2")
                         "both": Return "Broad (Split)" format when both exist (e.g., A*02:03 -> "A2 (A203)")
        
        Returns:
            Converted HLA type(s). Returns string for serological, list for molecular.
        """
        # Validate handle_broad parameter
        if handle_broad not in ["split", "broad", "both"]:
            raise ValueError(f"Invalid handle_broad value: {handle_broad}. Use 'split', 'broad', or 'both'.")
        
        if target_format is None:
            # Auto-detect and convert to opposite
            if self._is_serological(hla_type):
                return self._to_molecular(hla_type, expand_splits)
            else:
                return self._to_serology(hla_type, handle_broad)
        elif target_format == "s":
            return self._to_serology(hla_type, handle_broad)
        elif target_format == "m":
            return self._to_molecular(hla_type, expand_splits)
        else:
            raise ValueError(f"Unsupported target format: {target_format}. Use 's' for serological, 'm' for molecular, or None for auto-detect.")
    
    def _to_2field(self, molecular_allele):
        """Convert molecular allele to 2-field format (e.g., A*01:01:01:01 -> A*01:01)."""
        if '*' in molecular_allele and ':' in molecular_allele:
            parts = molecular_allele.split(':')
            if len(parts) >= 2:
                return f"{parts[0]}:{parts[1]}"
        
        # If it's not a molecular allele or doesn't have enough fields, return as-is
        return molecular_allele
    
    def _to_serology(self, hla_type, handle_broad="split"):
        """Convert molecular allele to serological equivalent."""
        if '*' in hla_type:
            # Validate that this is a real molecular allele first
            if not self._is_valid_molecular_allele(hla_type):
                raise HLAConversionError(f"Unrecognized molecular allele: {hla_type}")
            
            # Try exact match first
            sero = self.parser.get_molecular_to_serological(hla_type, return_broad=False)
            if not sero:
                # Try 2-field version
                two_field = self._to_2field(hla_type)
                sero = self.parser.get_molecular_to_serological(two_field, return_broad=False)
            
            if sero:
                return self._format_serological_result(sero, handle_broad)
            
            # If no serological mapping found
            raise HLAConversionError(f"No serological equivalent found for molecular allele: {hla_type}")
        
        # If already serological, validate and return as-is
        if self._is_serological(hla_type):
            if not self._is_valid_serological_allele(hla_type):
                raise HLAConversionError(f"Unrecognized serological allele: {hla_type}")
            return hla_type
        
        # Invalid format
        raise HLAConversionError(f"Invalid HLA format: {hla_type}")
    
    def _format_serological_result(self, sero, handle_broad):
        """Format serological result based on handle_broad setting."""
        if handle_broad == "split":
            return sero
        elif handle_broad == "broad":
            # Try to get broad antigen, return original if no broad exists
            broad = self.parser.get_broad_antigen(sero)
            return broad if broad else sero
        elif handle_broad == "both":
            # Get broad antigen if it exists
            broad = self.parser.get_broad_antigen(sero)
            if broad:
                return f"{broad} ({sero})"
            else:
                return sero
        else:
            return sero    
    
    def _to_molecular(self, hla_type, expand_splits=False):
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
    
    def _is_serological(self, hla_type):
        """Check if HLA type is in serological format."""
        # Serological types typically don't contain '*' and may contain numbers
        # Handle standard format (A1, B27) and C locus format (Cw14)
        return '*' not in hla_type and (re.match(r'^[A-Z]+\d+$', hla_type) or re.match(r'^Cw\d+$', hla_type))
    
    def _is_valid_molecular_allele(self, hla_type):
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
    
    def _is_valid_serological_allele(self, hla_type):
        """Check if serological allele is known."""
        # Ensure data is loaded
        self.parser._load_data()
        # Check if it's in our serological mapping
        return hla_type in self.parser._serological_mapping
    
    def classify_kir_ligand(self, hla_type, bead_annotation=None):
        """Classify HLA type for KIR ligand properties.
        
        Args:
            hla_type: HLA type (molecular or serological)
            bead_annotation: Optional bead annotation (e.g., "Bw6" from "B7,Bw6")
            
        Returns:
            Dictionary with KIR ligand classification results
        """
        if not self.kir_classifier:
            raise RuntimeError("KIR ligand classification not enabled. Initialize with enable_kir=True")
        
        # Ensure KIR data is loaded
        self.kir_classifier.load_data()
        
        # If molecular allele, classify directly
        if '*' in hla_type:
            return self.kir_classifier.classify_allele(hla_type, bead_annotation)
        
        # If serological, get molecular alleles and classify
        if self._is_serological(hla_type):
            try:
                molecular_alleles = self._to_molecular(hla_type)
                return self.kir_classifier.classify_serological(hla_type, molecular_alleles, bead_annotation)
            except HLAConversionError:
                # If conversion fails, return empty result
                return {
                    "is_kir_ligand": False,
                    "kir_ligand_type": None,
                    "kir_receptors": [],
                    "source": None
                }
        
        # Invalid format
        return {
            "is_kir_ligand": False,
            "kir_ligand_type": None,
            "kir_receptors": [],
            "source": None
        }


def convert(hla_type, version, target_format=None, data_dir="data", expand_splits=False, handle_broad="split"):
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
        handle_broad: How to handle broad/split relationships in molecular to serological conversion:
                     "split" (default): Return split antigen (e.g., A*02:03 -> "A203")
                     "broad": Return broad antigen when available (e.g., A*02:03 -> "A2")  
                     "both": Return "Broad (Split)" format when both exist (e.g., A*02:03 -> "A2 (A203)")
    
    Returns:
        Converted HLA type(s). Returns string for serological, list for molecular.
    """
    converter = HLAConverter(version, data_dir)
    return converter.convert(hla_type, target_format, expand_splits, handle_broad)