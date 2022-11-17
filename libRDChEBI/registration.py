from chembl_structure_pipeline.standardizer import parse_molblock
from rdkit.Chem import RegistrationHash


def get_registration_hash(molfile):
    mol = parse_molblock(molfile)
    layers = RegistrationHash.GetMolLayers(mol)
    r_hash = RegistrationHash.GetMolHash(layers)
    return r_hash
