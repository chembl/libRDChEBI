from chembl_structure_pipeline.standardizer import (
    parse_molblock,
    update_mol_valences,
)
from rdkit.Chem import Descriptors, Mol, Atom
from rdkit import Chem
from typing import Optional
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
        str: Molecular formula or None if unable to generate
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
        formula (str): Molecular formula
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


def _get_frag_formula(mol: Mol) -> str:
    """
    Generate molecular formula for a molecule fragment.
    Handles special cases like R groups, isotopes, and dummy atoms.

    Args:
        mol (rdkit.Chem.Mol): RDKit molecule object

    Returns:
        str: Molecular formula
    """
    atoms_dict = {}
    isotopes_dict = {}
    hs = 0

    # Count atoms and handle special cases
    for at in mol.GetAtoms():
        symbol = at.GetSymbol()
        isotope = at.GetIsotope()

        if isotope and at.GetAtomicNum() != 0 and symbol != "H":
            isotopes_dict.setdefault(symbol, {})
            isotopes_dict[symbol][isotope] = isotopes_dict[symbol].get(isotope, 0) + 1
        else:
            if atom_is_r_group(at):
                atoms_dict["R"] = atoms_dict.get("R", 0) + 1
            elif symbol == "H":
                if isotope == 2:
                    atoms_dict["D"] = atoms_dict.get("D", 0) + 1
                elif isotope == 3:
                    atoms_dict["T"] = atoms_dict.get("T", 0) + 1
                else:
                    hs += 1
            elif symbol == "*" and at.GetQueryType():
                atoms_dict[at.GetQueryType()] = atoms_dict.get(at.GetQueryType(), 0) + 1
            else:
                atoms_dict[symbol] = atoms_dict.get(symbol, 0) + 1

        hs += at.GetTotalNumHs(includeNeighbors=False)

    if hs > 0:
        atoms_dict["H"] = hs

    # Remove dummy atoms
    atoms_dict.pop("*", None)

    # Handle R groups
    r_part = (
        f"R{atoms_dict['R']}"
        if "R" in atoms_dict and atoms_dict["R"] > 1
        else "R"
        if "R" in atoms_dict
        else ""
    )
    atoms_dict.pop("R", None)

    # Build formula string following Hill notation
    elements = ["C", "H"] + sorted(set(atoms_dict.keys()) - {"C", "H"})
    formula = ""

    for elem in elements:
        if elem in atoms_dict:
            count = atoms_dict[elem]
            formula += elem if count == 1 else f"{elem}{count}"

        if elem in isotopes_dict:
            for iso, count in sorted(isotopes_dict[elem].items()):
                formula += f"[{iso}{elem}]" if count == 1 else f"[{iso}{elem}{count}]"

    return formula + r_part


def _get_small_molecule_formula(molfile: str) -> str:
    """
    Generate molecular formula for a small molecule.
    Handles multiple fragments and dummy atoms.

    Args:
        molfile (str): MOL file content as string

    Returns:
        str: Molecular formula fragments separated by dots
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
        str: Molecular formula with substance group labeling, or None if fails
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
