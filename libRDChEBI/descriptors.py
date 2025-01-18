from chembl_structure_pipeline.standardizer import (
    parse_molblock,
    update_mol_valences,
)
from rdkit.Chem import Descriptors, Mol, Atom
from rdkit import Chem
from typing import Dict, Optional, Union
import re

polymer_regex = re.compile(
    r"^M  STY.+(SRU)|(MON)|(COP)|(CRO)|(ANY)", flags=re.MULTILINE
)


def is_polymer(molfile: str) -> bool:
    """
    Check if the molecule is a polymer based on MOL file structure type flags.

    Args:
        molfile (str): MOL file content as string

    Returns:
        bool: True if molecule is a polymer, False otherwise
    """
    if polymer_regex.search(molfile):
        return True
    else:
        return False


def get_molformula(molfile: str) -> Optional[str]:
    """
    Get molecular formula for any molecule type (polymer or small molecule).

    Args:
        molfile (str): MOL file content as string

    Returns:
        str: Molecular formula following Hill notation, or None if unable to generate
    """
    if is_polymer(molfile):
        return _get_polymer_formula(molfile)
    else:
        return _get_small_molecule_formula(molfile)


def get_avg_mass(molfile: str) -> Optional[float]:
    """
    Calculate average molecular mass.

    Args:
        molfile (str): MOL file content as string

    Returns:
        float: Average molecular mass or None if calculation fails
    """
    avg_mass = None
    mol = parse_molblock(molfile)
    if mol:
        mol = update_mol_valences(mol)
        avg_mass = Descriptors.MolWt(mol)
    return avg_mass


def get_monoisotopic_mass(molfile: str) -> Optional[float]:
    """
    Calculate monoisotopic molecular mass.

    Args:
        molfile (str): MOL file content as string

    Returns:
        float: Monoisotopic mass or None if calculation fails
    """
    monoisotopic_mass = None
    mol = parse_molblock(molfile)
    if mol:
        mol = update_mol_valences(mol)
        monoisotopic_mass = Descriptors.ExactMolWt(mol)
    return monoisotopic_mass


def get_net_charge(molfile: str) -> int:
    """
    Calculate the net formal charge of a molecule.

    Args:
        molfile (str): MOL file content as string

    Returns:
        int: Sum of all formal charges in the molecule
    """
    mol = parse_molblock(molfile)
    charges = [atm.GetFormalCharge() for atm in mol.GetAtoms()]
    return sum(charges)


def get_mass_from_formula(formula: str, average: bool = True) -> Optional[float]:
    """
    Calculate molecular mass from a molecular formula string.

    Args:
        formula (str): Molecular formula in Hill notation
        average (bool): If True, calculate average mass; if False, calculate monoisotopic mass

    Returns:
        float: Calculated mass or None if formula is invalid
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
            func = periodic_table.GetAtomicWeight
        else:
            func = periodic_table.GetMostCommonIsotopeMass
        try:
            elem_mass = func(matches[idx])
        except RuntimeError as e:
            return None

        mass += elem_mass * mult
    return mass


def atom_is_r_group(at: Atom) -> bool:
    """
    Check if an atom represents an R group.
    Excludes actual elements starting with R (Ra, Rb, Re, Rf, Rg, Rh, Rn, Ru).

    Args:
        at (rdkit.Chem.Atom): RDKit atom object

    Returns:
        bool: True if atom is an R group, False otherwise
    """
    # we don't want to mess with Ra, Rb, Re, Rf, Rg, Rh, Rn, Ru
    # to make sure is an R grup (R, R#, R1, Rn... ) AtomicNum must be 0
    if at.GetSymbol()[0] == "R" and at.GetAtomicNum() == 0:
        return True
    else:
        return False


def has_r_group(molfile: str) -> bool:
    """
    Check if molecule contains any R groups.

    Args:
        molfile (str): MOL file content as string

    Returns:
        bool: True if molecule contains R groups, False otherwise
    """
    mol = parse_molblock(molfile)
    for at in mol.GetAtoms():
        if atom_is_r_group(at):
            return True
    return False


def no_r_group_and_alias(molfile: str) -> bool:
    """
    Check if molecule has R groups defined as Carbon atoms with R aliases.
    This is used to handle legacy ChEBI molecules with non-standard R group definitions.

    Args:
        molfile (str): MOL file content as string

    Returns:
        bool: True if molecule has R group aliases but no proper R groups, False otherwise
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


def has_dummy_atom(molfile: str) -> bool:
    """
    Check if molecule contains dummy atoms (*).

    Args:
        molfile (str): MOL file content as string

    Returns:
        bool: True if molecule contains dummy atoms, False otherwise
    """
    mol = parse_molblock(molfile)
    for at in mol.GetAtoms():
        if at.GetSymbol() == "*":
            return True
    return False


def _create_or_add_one(dictionary: Dict, key: str) -> Dict:
    """
    Increment value for key in dictionary, creating key if it doesn't exist.

    Args:
        dictionary (dict): Dictionary to modify
        key: Dictionary key

    Returns:
        dict: Modified dictionary
    """
    dictionary[key] = dictionary.get(key, 0) + 1
    return dictionary


def _create_or_add_one_nested(
    dictionary: Dict, key1: str, key2: Union[str, int]
) -> Dict:
    """
    Increment value in nested dictionary, creating keys if they don't exist.

    Args:
        dictionary (dict): Dictionary to modify
        key1: Outer dictionary key
        key2: Inner dictionary key

    Returns:
        dict: Modified dictionary
    """
    dictionary[key1] = dictionary.get(key1, {})
    dictionary[key1][key2] = dictionary[key1].get(key2, 0) + 1
    return dictionary


def _get_frag_formula(mol: Mol) -> str:
    """
    Generate molecular formula for a molecule fragment.
    Handles special cases like R groups, isotopes, and dummy atoms.

    Args:
        mol (rdkit.Chem.Mol): RDKit molecule object

    Returns:
        str: Molecular formula following Hill notation
    """
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


def _get_small_molecule_formula(molfile: str) -> str:
    """
    Generate molecular formula for a small molecule.
    Handles multiple fragments and dummy atoms.

    Args:
        molfile (str): MOL file content as string

    Returns:
        str: Molecular formula following Hill notation, fragments separated by dots
    """
    mol = parse_molblock(molfile)
    mol = update_mol_valences(mol)
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    formulas = [_get_frag_formula(frag) for frag in frags]
    # disconnected dummy atom woud generate '' as a formula.
    # don't want to concatenate that
    return ".".join(filter(None, formulas))


def _get_polymer_formula(molfile: str) -> Optional[str]:
    """
    Generate molecular formula for a polymer.
    Handles substance groups and remaining atoms outside groups.

    Args:
        molfile (str): MOL file content as string

    Returns:
        str: Molecular formula following Hill notation with substance group labeling, or None if fails
    """
    mol = parse_molblock(molfile)
    mol = update_mol_valences(mol)
    sgroups = Chem.GetMolSubstanceGroups(mol)
    formulas = []
    processed_atoms = set()

    # First pass - process all defined sgroups
    for sg in sgroups:
        if not sg.HasProp("TYPE"):
            continue

        sg_type = sg.GetProp("TYPE")
        if sg_type in ("SUP", "MUL"):
            continue

        sg_atoms = set(sg.GetAtoms())
        if not sg_atoms:
            continue

        # Create submolecule for this sgroup
        sg_mol = Chem.RWMol()
        atom_map = {}

        for at_idx in sg_atoms:
            if at_idx in processed_atoms:
                continue
            atom = mol.GetAtomWithIdx(at_idx)
            new_idx = sg_mol.AddAtom(atom)
            atom_map[at_idx] = new_idx
            processed_atoms.add(at_idx)

        # Add bonds between atoms in this sgroup
        for at_idx in sg_atoms:
            atom = mol.GetAtomWithIdx(at_idx)
            for bond in atom.GetBonds():
                begin_idx = bond.GetBeginAtomIdx()
                end_idx = bond.GetEndAtomIdx()
                if begin_idx in sg_atoms and end_idx in sg_atoms:
                    if begin_idx not in atom_map or end_idx not in atom_map:
                        continue
                    if not sg_mol.GetBondBetweenAtoms(
                        atom_map[begin_idx], atom_map[end_idx]
                    ):
                        sg_mol.AddBond(
                            atom_map[begin_idx], atom_map[end_idx], bond.GetBondType()
                        )

        # Get formula for this sgroup
        sg_formula = _get_frag_formula(sg_mol)
        if not sg_formula:
            continue

        # Add label if present
        label = ""
        if sg.HasProp("LABEL"):
            label = sg.GetProp("LABEL")

        formula = f"({sg_formula}){label}"
        formulas.append(formula)

    # Second pass - collect all remaining atoms into a single molecule
    remaining_mol = Chem.RWMol()
    atom_map = {}

    for i in range(mol.GetNumAtoms()):
        if i not in processed_atoms:
            atom = mol.GetAtomWithIdx(i)
            new_idx = remaining_mol.AddAtom(atom)
            atom_map[i] = new_idx

    for i in range(mol.GetNumBonds()):
        bond = mol.GetBondWithIdx(i)
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        if begin_idx in atom_map and end_idx in atom_map:
            remaining_mol.AddBond(
                atom_map[begin_idx], atom_map[end_idx], bond.GetBondType()
            )

    # Get formula for remaining atoms if any
    if remaining_mol.GetNumAtoms() > 0:
        remaining_formula = _get_frag_formula(remaining_mol)
        if remaining_formula:
            formulas.append(remaining_formula)

    if not formulas:
        return None

    return ".".join(formulas)
