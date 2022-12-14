from chembl_structure_pipeline.standardizer import parse_molblock
from rdkit.Chem import RegistrationHash


def get_registration_layers(molfile):
    mol = parse_molblock(molfile)
    layers = RegistrationHash.GetMolLayers(mol)
    return layers


def get_registration_hash(layers):
    r_hash = RegistrationHash.GetMolHash(layers)
    return r_hash
