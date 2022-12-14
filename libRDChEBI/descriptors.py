from chembl_structure_pipeline.standardizer import (
    parse_molblock,
    update_mol_valences,
    get_isotope_parent_mol,
)
from rdkit.Chem import Descriptors
from rdkit import Chem
import re


polymer_regex = re.compile(
    r"^M  STY.+(SRU)|(MON)|(COP)|(CRO)|(ANY)", flags=re.MULTILINE
)


def has_r_group(molfile):
    mol = parse_molblock(molfile)
    for at in mol.GetAtoms():
        if at.GetSymbol().startswith("R"):
            return True
    return False


def has_dummy_atom(molfile):
    mol = parse_molblock(molfile)
    for at in mol.GetAtoms():
        if at.GetSymbol() == '*':
            return True
    return False


def is_polymer(molfile):
    if polymer_regex.search(molfile):
        return True
    else:
        return False


def get_net_charge(molfile):
    mol = parse_molblock(molfile)
    charges = [atm.GetFormalCharge() for atm in mol.GetAtoms()]
    return sum(charges)


def _get_frag_formula(mol):
    atoms_dict = {}
    for at in mol.GetAtoms():
        if at.GetSymbol()[0] == "R":
            if atoms_dict.get("R"):
                atoms_dict["R"] += 1
            else:
                atoms_dict["R"] = 1
        else:
            if atoms_dict.get(at.GetSymbol()):
                atoms_dict[at.GetSymbol()] += 1
            else:
                atoms_dict[at.GetSymbol()] = 1
    hs = 0
    for at in mol.GetAtoms():
        if at.GetSymbol() == "H":
            hs += 1
        hs += at.GetTotalNumHs(includeNeighbors=False)
    if hs > 0:
        atoms_dict["H"] = hs

    # '*' represent fragments (attaches to something)
    # and do not appear in molecular formula in ChEBI
    if atoms_dict.get("*"):
        del atoms_dict["*"]

    # don't show the number of atoms if count is 1
    atom_str_counts = lambda x: f"{x}" if atoms_dict[x] == 1 else f"{x}{atoms_dict[x]}"

    # R represents a class of compounds (something attaches here)
    # it appears in the molecular formula
    rs = ""
    if atoms_dict.get("R"):
        rs = atom_str_counts("R")
        del atoms_dict["R"]

    # Hill order system: carbon, hydrogen, then all other elements in alphabetical order
    molecular_formula = ""
    for elem in ("C", "H"):
        if atoms_dict.get(elem):
            molecular_formula += atom_str_counts(elem)
            del atoms_dict[elem]
    for at in sorted(atoms_dict.keys()):
        molecular_formula += atom_str_counts(at)
    molecular_formula = molecular_formula + rs
    return molecular_formula


def get_small_molecule_formula(molfile):
    mol = parse_molblock(molfile)
    mol = update_mol_valences(mol)
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    formulas = [_get_frag_formula(frag) for frag in frags]
    return ".".join(formulas)


def get_avg_mass(molfile):
    avg_mass = None
    mol = parse_molblock(molfile)
    if mol:
        mol = update_mol_valences(mol)
        # When removing isotpe information following RDKit functions will
        #  - MolWt: Calc the average weight.
        #  - ExactMolWt: Calc the monoisotopic weight
        mol = get_isotope_parent_mol(mol)
        avg_mass = Descriptors.MolWt(mol)
    return avg_mass


def get_monoisotopic_mass(molfile):
    monoisotopic_mass = None
    mol = parse_molblock(molfile)
    if mol:
        mol = update_mol_valences(mol)
        mol = get_isotope_parent_mol(mol)
        monoisotopic_mass = Descriptors.ExactMolWt(mol)
    return monoisotopic_mass


def get_polymer_formula(molfile):
    mol = parse_molblock(molfile)
    mol = update_mol_valences(mol)
    rwmol = Chem.RWMol(mol)
    formulas = []
    atoms_in_sgroups = []
    for sg in Chem.GetMolSubstanceGroups(rwmol):
        sub_mol = Chem.RWMol()
        # we only need the atoms (with their valences) in SGroups for the formula
        for atm in sg.GetAtoms():
            atom = rwmol.GetAtomWithIdx(atm)
            sub_mol.AddAtom(atom)
            atoms_in_sgroups.append(atm)

        formula = _get_frag_formula(sub_mol)
        if sg.HasProp("LABEL"):
            label = sg.GetProp("LABEL")
        else:
            label = ""
        formula = f"({formula}){label}"
        formulas.append(formula)

    # calc formula for the rest of atoms
    rwmol.BeginBatchEdit()
    for atm in atoms_in_sgroups:
        rwmol.RemoveAtom(atm)
    rwmol.CommitBatchEdit()
    rest_formula = _get_frag_formula(rwmol)

    if rest_formula:
        formulas.append(rest_formula)
    return ".".join(formulas)


def get_polymer_mass(molfile, avg=True):
    if avg:
        func = Descriptors.MolWt
    else:
        func = Descriptors.ExactMolWt
    mol = parse_molblock(molfile)
    mol = update_mol_valences(mol)
    rwmol = Chem.RWMol(mol)
    masses = []
    atoms_in_sgroups = []
    for sg in Chem.GetMolSubstanceGroups(rwmol):
        sub_mol = Chem.RWMol()
        for atm in sg.GetAtoms():
            atom = rwmol.GetAtomWithIdx(atm)
            sub_mol.AddAtom(atom)
            atoms_in_sgroups.append(atm)

        mass = round(func(sub_mol), 5)
        if sg.HasProp("LABEL"):
            label = sg.GetProp("LABEL")
        else:
            label = ""
        mass = f"({mass}){label}"
        masses.append(mass)

    # calc the mass for the rest of atoms
    rwmol.BeginBatchEdit()
    for atm in atoms_in_sgroups:
        rwmol.RemoveAtom(atm)
    rwmol.CommitBatchEdit()
    rest_mass = round(func(rwmol), 5)
    if rest_mass > 0.0:  # potential remaining '*' have mass 0.0
        masses.append(str(rest_mass))
    return "+".join(masses)


def get_mass_from_formula(formula, average=True):
    """
    average=True: avg mass
    average=False: monoisotopic mass
    """
    periodic_table = Chem.GetPeriodicTable()
    matches = re.findall("[A-Z][a-z]?|[0-9]+", formula)
    mass = 0
    for idx in range(len(matches)):
        if matches[idx].isnumeric():
            continue
        mult = (
            int(matches[idx + 1])
            if len(matches) > idx + 1 and matches[idx + 1].isnumeric()
            else 1
        )
        if average:
            elem_mass = periodic_table.GetAtomicWeight(matches[idx])
        else:
            elem_mass = periodic_table.GetMostCommonIsotopeMass(matches[idx])
        mass += elem_mass * mult
    return round(mass, 5)
