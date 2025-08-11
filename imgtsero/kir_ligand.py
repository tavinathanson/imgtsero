"""KIR ligand classification functionality for HLA alleles."""

import os
import json
import urllib.request
import urllib.parse
import urllib.error
from collections import defaultdict


class KIRLigandClassifier:
    """Classifier for HLA KIR ligand types."""
    
    def __init__(self, version, data_dir="data"):
        """Initialize KIR ligand classifier.
        
        Args:
            version: IPD-IMGT/HLA version (e.g., "3.61.0" or 3610)
            data_dir: Directory to store cached KIR ligand data
        """
        # Convert version format if needed
        self.version = self._normalize_version(version)
        self.data_dir = data_dir
        self._kir_ligand_map = {}  # Maps allele name -> KIR ligand type
        self._loaded = False
        
        # KIR receptor mappings
        self.kir_receptors = {
            "Bw4": ["KIR3DL1"],
            "C1": ["KIR2DL2", "KIR2DL3"],
            "C2": ["KIR2DL1"]
        }
        
        # Ensure data directory exists
        if not os.path.exists(data_dir):
            os.makedirs(data_dir)
    
    def _normalize_version(self, version):
        """Convert version format between 3610 and 3.61.0 formats."""
        version_str = str(version)
        
        # If it's already in X.XX.X format, return as-is
        if '.' in version_str and version_str.count('.') == 2:
            return version_str
        
        # Convert from 3610 to 3.61.0 format
        if len(version_str) == 4 and version_str.isdigit():
            major = version_str[0]
            minor = version_str[1:3].lstrip('0') or '0'
            return f"{major}.{minor}.0"
        
        # Try to parse other formats
        try:
            # Remove any trailing .0 if present
            if version_str.endswith('.0'):
                version_str = version_str[:-2]
            
            parts = version_str.split('.')
            if len(parts) == 2:
                return f"{parts[0]}.{parts[1]}.0"
            elif len(parts) == 1 and len(parts[0]) <= 2:
                # Single number like "3" -> "3.0.0"
                return f"{parts[0]}.0.0"
        except:
            pass
        
        # Return as-is if we can't parse it
        return version_str
    
    def _fetch_kir_data_from_api(self):
        """Fetch KIR ligand data from IPD-IMGT/HLA API.
        
        Returns:
            Dictionary mapping allele names to KIR ligand types
        """
        base_url = "https://www.ebi.ac.uk/cgi-bin/ipd/api/allele"
        
        # Build query parameters
        params = {
            'fields': 'name,locus,matching.kir_ligand,release_version',
            'query': f'and(or(eq(locus,"B*"),eq(locus,"C*")), eq(release_version,"{self.version}"))',
            'limit': '100000',
            'format': 'json'
        }
        
        # Construct URL with encoded parameters
        query_string = urllib.parse.urlencode(params)
        url = f"{base_url}?{query_string}"
        
        try:
            # Fetch data from API
            with urllib.request.urlopen(url) as response:
                data = json.loads(response.read().decode('utf-8'))
                
            if not data.get('data'):
                raise RuntimeError(f"No KIR ligand data found for version {self.version}")
            
            # Process API response
            kir_map = {}
            for entry in data['data']:
                allele_name = entry.get('name', '')
                kir_ligand = entry.get('matching.kir_ligand')
                if allele_name:
                    kir_map[allele_name] = kir_ligand
                    
            return kir_map
            
        except urllib.error.HTTPError as e:
            if e.code == 404:
                raise RuntimeError(f"KIR ligand data not available for version {self.version}")
            else:
                raise RuntimeError(f"Error fetching KIR ligand data: HTTP {e.code}")
        except urllib.error.URLError as e:
            raise RuntimeError(f"Error connecting to IPD-IMGT/HLA API: {e}")
        except json.JSONDecodeError as e:
            raise RuntimeError(f"Error parsing API response: {e}")
    
    def _get_cache_filename(self):
        """Get the cache filename for this version."""
        return os.path.join(self.data_dir, f"kir_ligand_{self.version}.json")
    
    def _load_from_cache(self):
        """Try to load KIR ligand data from cache.
        
        Returns:
            True if loaded successfully, False otherwise
        """
        cache_file = self._get_cache_filename()
        if os.path.exists(cache_file):
            try:
                with open(cache_file, 'r') as f:
                    self._kir_ligand_map = json.load(f)
                return True
            except (json.JSONDecodeError, IOError):
                return False
        return False
    
    def _save_to_cache(self):
        """Save KIR ligand data to cache."""
        cache_file = self._get_cache_filename()
        try:
            with open(cache_file, 'w') as f:
                json.dump(self._kir_ligand_map, f, indent=2)
        except IOError as e:
            print(f"Warning: Could not save KIR ligand cache: {e}")
    
    def load_data(self, force_download=False):
        """Load KIR ligand data, fetching from API if necessary.
        
        Args:
            force_download: Force download from API even if cache exists
        """
        if self._loaded and not force_download:
            return
        
        # Try to load from cache first
        if not force_download and self._load_from_cache():
            print(f"Loaded KIR ligand data from cache for version {self.version}")
            self._loaded = True
            return
        
        # Fetch from API
        print(f"Fetching KIR ligand data from IPD-IMGT/HLA API for version {self.version}...")
        self._kir_ligand_map = self._fetch_kir_data_from_api()
        
        # Save to cache
        self._save_to_cache()
        print(f"Successfully fetched and cached KIR ligand data for version {self.version}")
        
        self._loaded = True
    
    def get_kir_ligand(self, allele_name):
        """Get KIR ligand type for an allele.
        
        Args:
            allele_name: HLA allele name (e.g., "B*07:02:01")
            
        Returns:
            KIR ligand type (e.g., "Bw4", "Bw6", "C1", "C2") or None
        """
        if not self._loaded:
            self.load_data()
        
        # Direct lookup
        if allele_name in self._kir_ligand_map:
            return self._kir_ligand_map[allele_name]
        
        # Try progressively shorter forms (remove fields from right)
        parts = allele_name.split(':')
        for i in range(len(parts) - 1, 0, -1):
            shorter_name = ':'.join(parts[:i])
            if shorter_name in self._kir_ligand_map:
                return self._kir_ligand_map[shorter_name]
        
        return None
    
    def classify_allele(self, allele_name, bead_annotation=None):
        """Classify an allele for KIR ligand properties.
        
        Args:
            allele_name: HLA allele name (e.g., "B*07:02:01")
            bead_annotation: Optional bead annotation for validation (e.g., "Bw6" from "B7,Bw6")
            
        Returns:
            Dictionary with classification results:
                - is_kir_ligand: Boolean indicating if it's a KIR ligand
                - kir_ligand_type: The ligand type (Bw4, C1, C2) or None
                - kir_receptors: List of KIR receptors that bind this ligand
                - source: Always "api" since we use API as source of truth
        """
        # Always get classification from API data
        kir_ligand = self.get_kir_ligand(allele_name)
        
        result = {
            "is_kir_ligand": False,
            "kir_ligand_type": kir_ligand,
            "kir_receptors": [],
            "source": "api" if kir_ligand else None
        }
        
        if kir_ligand:
            result["is_kir_ligand"] = kir_ligand in ["Bw4", "C1", "C2"]
            if result["is_kir_ligand"]:
                result["kir_receptors"] = self.kir_receptors.get(kir_ligand, [])
        
        # Validate against bead annotation if provided
        if bead_annotation and bead_annotation in ["Bw4", "Bw6"]:
            if kir_ligand and kir_ligand != bead_annotation:
                raise ValueError(
                    f"Bead annotation '{bead_annotation}' conflicts with API data '{kir_ligand}' "
                    f"for allele {allele_name}"
                )
        
        return result
    
    def classify_serological(self, serological_name, molecular_alleles, bead_annotation=None):
        """Classify a serological antigen based on its molecular alleles.
        
        Args:
            serological_name: Serological antigen name (e.g., "B7", "Cw7")
            molecular_alleles: List of molecular alleles for this serological type
            bead_annotation: Optional bead annotation for validation (e.g., "Bw6")
            
        Returns:
            Dictionary with classification results (same format as classify_allele)
        """
        # Always check molecular alleles from API
        kir_ligand_types = set()
        for allele in molecular_alleles:
            kir_ligand = self.get_kir_ligand(allele)
            if kir_ligand:
                kir_ligand_types.add(kir_ligand)
        
        # Determine consensus KIR ligand type
        if len(kir_ligand_types) == 1:
            kir_ligand = kir_ligand_types.pop()
            result = {
                "is_kir_ligand": kir_ligand in ["Bw4", "C1", "C2"],
                "kir_ligand_type": kir_ligand,
                "kir_receptors": self.kir_receptors.get(kir_ligand, []) if kir_ligand in ["Bw4", "C1", "C2"] else [],
                "source": "api"
            }
        elif len(kir_ligand_types) > 1:
            # Mixed types - this shouldn't happen for a proper serological type
            result = {
                "is_kir_ligand": False,
                "kir_ligand_type": f"Mixed: {', '.join(sorted(kir_ligand_types))}",
                "kir_receptors": [],
                "source": "api"
            }
        else:
            # No KIR ligand data found
            result = {
                "is_kir_ligand": False,
                "kir_ligand_type": None,
                "kir_receptors": [],
                "source": None
            }
        
        # Validate against bead annotation if provided
        if bead_annotation and bead_annotation in ["Bw4", "Bw6"]:
            if result["kir_ligand_type"] and result["kir_ligand_type"] != bead_annotation:
                # Handle mixed types specially
                if not result["kir_ligand_type"].startswith("Mixed:"):
                    raise ValueError(
                        f"Bead annotation '{bead_annotation}' conflicts with API data '{result['kir_ligand_type']}' "
                        f"for serological antigen {serological_name}"
                    )
        
        return result
    
    def get_all_kir_ligands(self):
        """Get all alleles grouped by KIR ligand type.
        
        Returns:
            Dictionary mapping KIR ligand types to lists of alleles
        """
        if not self._loaded:
            self.load_data()
        
        grouped = defaultdict(list)
        for allele, kir_ligand in self._kir_ligand_map.items():
            if kir_ligand:
                grouped[kir_ligand].append(allele)
        
        return dict(grouped)
    
    @staticmethod
    def parse_bead_annotation(bead_name):
        """Parse SAB bead name to extract HLA antigen and KIR annotation.
        
        Args:
            bead_name: SAB bead name (e.g., "B27,Bw4" or "B7,Bw6" or "Cw7")
            
        Returns:
            Tuple of (HLA antigen, KIR annotation or None)
        """
        if ',' in bead_name:
            parts = bead_name.split(',')
            if len(parts) == 2:
                hla_antigen = parts[0].strip()
                annotation = parts[1].strip()
                # Only return annotation if it's Bw4 or Bw6
                if annotation in ["Bw4", "Bw6"]:
                    return hla_antigen, annotation
                return hla_antigen, None
        return bead_name.strip(), None