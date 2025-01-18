from chembl_structure_pipeline.standardizer import parse_molblock
from rdkit.Chem import RegistrationHash
from typing import Dict, Any


def get_registration_layers(molfile: str) -> Dict[str, Any]:
    """
    Generate registration layers for a molecule from its molfile.

    Args:
        molfile (str): The molecule structure in molfile format

    Returns:
        Dict[str, Any]: Dictionary containing the registration layers
    """
    mol = parse_molblock(molfile)
    layers = RegistrationHash.GetMolLayers(mol)
    return layers


def get_registration_hash(layers: Dict[str, Any]) -> str:
    """
    Generate a registration hash from molecule layers.

    Args:
        layers (Dict[str, Any]): Dictionary containing the registration layers

    Returns:
        str: The registration hash string
    """
    r_hash = RegistrationHash.GetMolHash(layers)
    return r_hash
