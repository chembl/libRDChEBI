from chembl_structure_pipeline.standardizer import (
    parse_molblock,
    update_mol_valences,
)
from rdkit.Chem import Descriptors
from rdkit import Chem
import re


polymer_regex = re.compile(
    r"^M  STY.+(SRU)|(MON)|(COP)|(CRO)|(ANY)", flags=re.MULTILINE
)


def atom_is_r_group(at):
    # we don't want to mess with Ra, Rb, Re, Rf, Rg, Rh, Rn, Ru
    # to make sure is an R grup (R, R#, R1, Rn... ) AtomicNum must be 0
    if at.GetSymbol()[0] == "R" and at.GetAtomicNum() == 0:
        return True
    else:
        return False


def has_r_group(molfile):
    mol = parse_molblock(molfile)
    for at in mol.GetAtoms():
        if atom_is_r_group(at):
            return True
    return False


def no_r_group_and_alias(molfile):
    """
    Some molecules in old ChEBI have R groups defined as Carbons with aliases
    this function finds them
    """
    mol = parse_molblock(molfile)
    r_group = False
    alias = False
    for at in mol.GetAtoms():
        if atom_is_r_group(at):
            r_group = True
        if "molFileAlias" in at.GetPropNames() and at.GetSymbol() == "C":
            if at.GetProp("molFileAlias")[0] == "R":
                alias = True
    if not r_group and alias:
        return True
    else:
        return False


def has_dummy_atom(molfile):
    mol = parse_molblock(molfile)
    for at in mol.GetAtoms():
        if at.GetSymbol() == "*":
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


def _create_or_add_one(dictionary, key):
    dictionary[key] = dictionary.get(key, 0) + 1
    return dictionary


def _create_or_add_one_nested(dictionary, key1, key2):
    dictionary[key1] = dictionary.get(key1, {})
    dictionary[key1][key2] = dictionary[key1].get(key2, 0) + 1
    return dictionary


def _get_frag_formula(mol):
    atoms_dict = {}
    isotopes_dict = {}
    hs = 0
    for at in mol.GetAtoms():
        # R groups can show isotopes(!) if numbered (e.g. R1, R2...)
        # so we want to exclude isotopes with AtomicNum=0 here
        # Deuterium and Tritium are a special case so skipping H here
        if at.GetIsotope() and at.GetAtomicNum() != 0 and at.GetSymbol() != "H":
            isotopes_dict = _create_or_add_one_nested(
                isotopes_dict, at.GetSymbol(), at.GetIsotope()
            )
        else:
            if atom_is_r_group(at):
                atoms_dict = _create_or_add_one(atoms_dict, "R")
            elif at.GetSymbol() == "H":
                if at.GetIsotope() == 2:  # Deuterium
                    atoms_dict = _create_or_add_one(atoms_dict, "D")
                elif at.GetIsotope() == 3:  # Tritium
                    atoms_dict = _create_or_add_one(atoms_dict, "T")
                else:
                    hs += 1
            # capture Reaxys generics
            elif at.GetSymbol() == "*" and at.GetQueryType():
                atoms_dict = _create_or_add_one(atoms_dict, at.GetQueryType())
            else:
                atoms_dict = _create_or_add_one(atoms_dict, at.GetSymbol())
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

    # Hill System Order:
    # Carbon first, Hydrogen second, and all remaining elements, including Deuterium and Tritium, in alphabetical order
    pse = Chem.GetPeriodicTable()
    els = list(
        filter(
            lambda x: x not in ("H", "C"),
            [pse.GetElementSymbol(i) for i in range(1, 119)],
        )
    )
    elements_list = ["C", "H"] + sorted(els + ["D", "T"]) + ["A", "M", "X"]

    molecular_formula = ""
    for elem in elements_list:
        if atoms_dict.get(elem):
            molecular_formula += atom_str_counts(elem)
        if isotopes_dict.get(elem):
            for iso in sorted(isotopes_dict[elem].keys()):
                if isotopes_dict[elem][iso] == 1:
                    molecular_formula += f"[{iso}{elem}]"
                else:
                    molecular_formula += f"[{iso}{elem}{isotopes_dict[elem][iso]}]"
    molecular_formula = molecular_formula + rs
    return molecular_formula


def get_small_molecule_formula(molfile):
    mol = parse_molblock(molfile)
    mol = update_mol_valences(mol)
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    formulas = [_get_frag_formula(frag) for frag in frags]
    # disconnected dummy atom woud generate '' as a formula.
    # don't want to concatenate that
    return ".".join(filter(None, formulas))


def get_avg_mass(molfile):
    avg_mass = None
    mol = parse_molblock(molfile)
    if mol:
        mol = update_mol_valences(mol)
        avg_mass = Descriptors.MolWt(mol)
    return avg_mass


def get_monoisotopic_mass(molfile):
    monoisotopic_mass = None
    mol = parse_molblock(molfile)
    if mol:
        mol = update_mol_valences(mol)
        monoisotopic_mass = Descriptors.ExactMolWt(mol)
    return monoisotopic_mass


def get_polymer_formula(molfile):
    mol = parse_molblock(molfile)
    mol = update_mol_valences(mol)
    rwmol = Chem.RWMol(mol)
    formulas = []
    atoms_to_remove = []

    for sg in Chem.GetMolSubstanceGroups(rwmol):
        sg_props = sg.GetPropsAsDict()
        if sg_props["TYPE"] in ("SUP", "MUL"):
            continue

        atoms_in_sgroup = [at_idx for at_idx in sg.GetAtoms()]
        atoms_to_remove += atoms_in_sgroup
        sg_sub_mol = Chem.RWMol()

        for at_idx in atoms_in_sgroup:
            atom = rwmol.GetAtomWithIdx(at_idx)
            sg_sub_mol.AddAtom(atom)

        sg_formula = _get_frag_formula(sg_sub_mol)
        conn_atoms = get_conn_atoms(rwmol, atoms_in_sgroup[0])
        remain_formula = ""

        if len(conn_atoms) > len(atoms_in_sgroup):
            sub_mol_conn = Chem.RWMol()
            for at_idx in set(conn_atoms) - set(atoms_in_sgroup):
                atom = rwmol.GetAtomWithIdx(at_idx)
                for ssg in Chem.GetMolSubstanceGroups(rwmol):
                    index = int(sg_props["index"])
                    if int(ssg.GetProp("index")) == index:
                        continue
                    if at_idx not in ssg.GetAtoms():
                        sub_mol_conn.AddAtom(atom)
                atoms_to_remove.append(at_idx)
            remain_formula = _get_frag_formula(sub_mol_conn)

        label = sg.GetProp("LABEL") if sg.HasProp("LABEL") else ""

        if remain_formula:
            formula = f"({sg_formula}){label}{remain_formula}"
        else:
            formula = f"({sg_formula}){label}"
        formulas.append(formula)

    rwmol.BeginBatchEdit()
    for atm in set(atoms_to_remove):
        rwmol.RemoveAtom(atm)
    rwmol.CommitBatchEdit()

    frags = Chem.GetMolFrags(rwmol, asMols=True, sanitizeFrags=False)
    remain_formulas = [_get_frag_formula(frag) for frag in frags]
    if remain_formulas:
        formulas += remain_formulas
    return ".".join(formulas)


def get_conn_atoms(mol, atom_idx):
    connected_atoms = set()
    queue = [atom_idx]
    visited = set()

    while queue:
        current_atom_idx = queue.pop(0)
        if current_atom_idx not in visited:
            visited.add(current_atom_idx)
            connected_atoms.add(current_atom_idx)
            neighbors = mol.GetAtomWithIdx(current_atom_idx).GetNeighbors()
            for neighbor in neighbors:
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in visited:
                    queue.append(neighbor_idx)
    return connected_atoms


def get_mass_from_formula(formula, average=True):
    """
    average=True: avg mass
    average=False: monoisotopic mass
    """
    periodic_table = Chem.GetPeriodicTable()
    matches = re.findall("[A-Z][a-z]?|[0-9]+", formula)
    mass = 0
    for idx in range(len(matches)):

        # skip R groups
        if matches[idx] == "R":
            continue

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
    return mass
