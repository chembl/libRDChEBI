from chembl_structure_pipeline.standardizer import parse_molblock, update_mol_valences
from rdkit.Chem import Descriptors
from rdkit import Chem
import ctfile
import io


def drop_isotopes_info(mol):
    """
    When removing isotpe information following RDKit functions will:
     - MolWt: Calc the average weight.
     - ExactMolWt: Calc the monoisotopic weight
    MolWt takes average weight of each atom only if no isotope information is given.
    ExactMolWt takes most abundant isotope weight for each atom only if no isotope information is given.
    """
    for atom in mol.GetAtoms():
        atom.SetIsotope(0)
    return mol


def get_net_charge(molfile):
    mol = parse_molblock(molfile)
    charges = [atm.GetFormalCharge() for atm in mol.GetAtoms()]
    return sum(charges)


def get_small_mol_formula(mol, ctf):
    atoms_dict = {}
    for rd_at in mol.GetAtoms():
        idx = rd_at.GetIdx()
        at = ctf.atoms[idx]
        if at.atom_symbol[0] == "R" or rd_at.GetSymbol().startswith("R"):
            if atoms_dict.get("R"):
                atoms_dict["R"] += 1
            else:
                atoms_dict["R"] = 1
        else:
            if atoms_dict.get(rd_at.GetSymbol()):
                atoms_dict[rd_at.GetSymbol()] += 1
            else:
                atoms_dict[rd_at.GetSymbol()] = 1

    hs = 0
    for at in mol.GetAtoms():
        if at.GetSymbol() == "H":
            hs += 1
        hs += at.GetTotalNumHs(includeNeighbors=False)
    if hs > 0:
        atoms_dict["H"] = hs

    # '*' represent groups and do not appear in molecular formula in ChEBI
    if atoms_dict.get("*"):
        del atoms_dict["*"]

    # don't show the number of atoms if count is 1
    atom_str_counts = lambda x: f"{x}" if atoms_dict[x] == 1 else f"{x}{atoms_dict[x]}"

    # R, in ChEBI, represents a class of compounds so it should appear in the molecular formula
    rr = ""
    if atoms_dict.get("R"):
        rr = atom_str_counts("R")
        del atoms_dict["R"]

    # Hill order system: carbon, hydrogen, then all other elements in alphabetical order
    molecular_formula = ""
    for elem in ("C", "H"):
        if atoms_dict.get(elem):
            molecular_formula += atom_str_counts(elem)
            del atoms_dict[elem]
    for k in sorted(atoms_dict.keys()):
        molecular_formula += atom_str_counts(k)
    molecular_formula = molecular_formula + rr
    return molecular_formula


def get_molecular_formula(molfile):
    mol = parse_molblock(molfile)
    mol = update_mol_valences(mol)
    f = io.StringIO(molfile)
    ctf = ctfile.load(f)

    rwmol = Chem.RWMol(mol)
    sgroups = Chem.GetMolSubstanceGroups(rwmol)

    formulas = []
    atoms_in_sgroups = []
    for sg in sgroups:
        sub_mol = Chem.RWMol()
        # we only need the atoms (with calculated valences) in SGroups for the formula
        for atm in sg.GetAtoms():
            atom = rwmol.GetAtomWithIdx(atm)
            sub_mol.AddAtom(atom)
            atoms_in_sgroups.append(atm)

        formula = ""
        formula = get_small_mol_formula(sub_mol, ctf)

        if sg.HasProp("LABEL"):
            label = sg.GetProp("LABEL")
        else:
            label = ""
        formula = f"({formula}){label}"
        formula = formula.replace("*", "")
        formulas.append(formula)

    # calc formula for the rest of atoms
    rwmol.BeginBatchEdit()
    for atm in atoms_in_sgroups:
        rwmol.RemoveAtom(atm)
    rwmol.CommitBatchEdit()
    rest_formula = get_small_mol_formula(rwmol, ctf)

    if rest_formula:
        formulas.append(rest_formula)
    mol_formula = ".".join(formulas)
    return mol_formula


def get_avg_mass(molfile):
    mol = parse_molblock(molfile)
    mol = update_mol_valences(mol)
    mol = drop_isotopes_info(mol)
    avg_mass = Descriptors.MolWt(mol)
    return avg_mass


def get_monoisotopic_mass(molfile):
    mol = parse_molblock(molfile)
    mol = update_mol_valences(mol)
    mol = drop_isotopes_info(mol)
    monoisotopic_mass = Descriptors.ExactMolWt(mol)
    return monoisotopic_mass


def get_inchi_and_key(molfile):
    inchi = Chem.MolBlockToInchi(molfile)
    inchi_key = Chem.InchiToInchiKey(inchi)
    return inchi, inchi_key


def get_smiles(molfile):
    # we do apply RDKit chemistry model for SMILES generation
    mol = Chem.MolFromMolBlock(molfile)
    smiles = Chem.MolToSmiles(mol)
    return smiles
