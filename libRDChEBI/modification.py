from chembl_structure_pipeline.standardizer import parse_molblock, update_mol_valences
from rdkit import Chem


def transform_alias_to_r(molfile):
    """
    Some molecules in old ChEBI have R groups defined as Carbons with aliases
    this function fixes them
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


def remove_hs(molfile):
    """Bespoke remove Hs function for MetaboLights team"""
    mol = parse_molblock(molfile)
    Chem.FastFindRings(mol)
    mol = update_mol_valences(mol)
    indices = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1 and not atom.GetIsotope():
            bnd = atom.GetBonds()[0]
            if not (
                bnd.GetBondDir() in (Chem.BondDir.BEGINWEDGE, Chem.BondDir.BEGINDASH)
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
