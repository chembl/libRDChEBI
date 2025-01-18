from chembl_structure_pipeline.standardizer import parse_molblock, update_mol_valences
from rdkit import Chem


def transform_alias_to_r(molfile: str) -> str:
    """
    Convert Carbon atoms with R-group aliases to proper R-group representations.

    Some molecules in old ChEBI have R groups defined as Carbons with aliases.
    This function converts them to proper R-group representations.

    Args:
        molfile (str): The molecule structure in molfile format

    Returns:
        str: Modified molecule structure in molfile format
    """
    mol = parse_molblock(molfile)
    for at in mol.GetAtoms():
        if "molFileAlias" in at.GetPropNames() and at.GetSymbol() == "C":
            alias = at.GetProp("molFileAlias")
            if alias.startswith("R"):
                at.SetAtomicNum(0)
                at.SetProp("dummyLabel", alias)
                at.SetProp("molFileAlias", "")
    return Chem.MolToMolBlock(mol)


def remove_hs(molfile: str) -> str:
    """
    Remove hydrogen atoms from a molecule structure.

    Bespoke remove Hs function for MetaboLights team. Preserves stereochemistry-relevant
    hydrogen atoms.

    Args:
        molfile (str): The molecule structure in molfile format

    Returns:
        str: Modified molecule structure in molfile format with hydrogens removed
    """
    mol = parse_molblock(molfile)
    Chem.FastFindRings(mol)
    mol = update_mol_valences(mol)
    indices = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1 and not atom.GetIsotope():
            bnd = atom.GetBonds()[0]
            if (
                bnd.GetBondDir()
                not in (Chem.BondDir.BEGINWEDGE, Chem.BondDir.BEGINDASH)
            ) and not (
                bnd.HasProp("_MolFileBondStereo")
                and bnd.GetUnsignedProp("_MolFileBondStereo") in (1, 6)
            ):
                indices.append(atom.GetIdx())
    mol = Chem.RWMol(mol)
    for index in sorted(indices, reverse=True):
        mol.RemoveAtom(index)
    props = molfile.split("M  END")[1].strip()
    props = props if len(props) > 1 else None
    out_molfile = Chem.MolToMolBlock(mol)
    if props:
        out_molfile += props
    return out_molfile
