from chembl_structure_pipeline.standardizer import parse_molblock, update_mol_valences
from rdkit import Chem


def get_smiles(molfile):
    mol = parse_molblock(molfile)
    # Minimal sanitization to generate SMILES
    if mol:
        mol.ClearComputedProps()
        mol = update_mol_valences(mol)
        Chem.SanitizeMol(
            mol,
            sanitizeOps=Chem.SANITIZE_SYMMRINGS
            | Chem.SANITIZE_FINDRADICALS
            | Chem.SANITIZE_SETAROMATICITY
            | Chem.SANITIZE_ADJUSTHS,
        )
    return Chem.MolToSmiles(mol) if mol else None
