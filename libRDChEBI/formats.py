from chembl_structure_pipeline.standardizer import parse_molblock, update_mol_valences
from rdkit import Chem


def get_smiles(molfile):
    mol = parse_molblock(molfile)
    mol = update_mol_valences(mol)
    Chem.SanitizeMol(
        mol,
        sanitizeOps=Chem.SanitizeFlags.SANITIZE_FINDRADICALS
        | Chem.SanitizeFlags.SANITIZE_KEKULIZE
        | Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
        | Chem.SanitizeFlags.SANITIZE_SETCONJUGATION
        | Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION
        | Chem.SanitizeFlags.SANITIZE_SYMMRINGS,
        catchErrors=True,
    )
    return Chem.MolToSmiles(mol)
