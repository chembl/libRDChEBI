from .registration import get_registration_hash
from .depiction import depict
from .descriptors import (
    get_smiles,
    get_avg_mass,
    get_molecular_formula,
    get_inchi_and_key,
    get_monoisotopic_mass,
    get_net_charge,
)
import rdkit

rdkversion = rdkit.__version__.split(".")
if rdkversion < ["2022", "09", "2"]:
    raise ValueError("need an RDKit version >= 2022.09.2")
