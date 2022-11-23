from chembl_structure_pipeline.standardizer import parse_molblock, update_mol_valences
from rdkit.Chem import Descriptors
from rdkit import Chem


def has_r_group(molfile):
    mol = parse_molblock(molfile)
    for at in mol.GetAtoms():
        if at.GetSymbol().startswith("R"):
            return True
    return False


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


def get_frag_formula(mol):
    atoms_dict = {}
    for at in mol.GetAtoms():
        if at.GetSymbol().startswith("R"):
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

    # '*' represent fragments and do not appear in molecular formula in ChEBI
    if atoms_dict.get("*"):
        del atoms_dict["*"]

    # don't show the number of atoms if count is 1
    atom_str_counts = lambda x: f"{x}" if atoms_dict[x] == 1 else f"{x}{atoms_dict[x]}"

    # R, in ChEBI, represents a class of compounds so it should appear in the molecular formula
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


def get_small_molecule_formula(mol):
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    formulas = [get_frag_formula(frag) for frag in frags]
    return ".".join(formulas)


def get_molecular_formula(molfile):
    mol = parse_molblock(molfile)
    mol = update_mol_valences(mol)
    if Chem.GetMolSubstanceGroups(mol):
        return get_polymer_formula(mol)
    else:
        return get_small_molecule_formula(mol)


def get_avg_mass(molfile):
    mol = parse_molblock(molfile)
    mol = update_mol_valences(mol)
    mol = drop_isotopes_info(mol)
    if Chem.GetMolSubstanceGroups(mol):
        return get_polymer_mass(mol, Descriptors.MolWt)
    else:
        return Descriptors.MolWt(mol)


def get_monoisotopic_mass(molfile):
    mol = parse_molblock(molfile)
    mol = update_mol_valences(mol)
    mol = drop_isotopes_info(mol)
    if Chem.GetMolSubstanceGroups(mol):
        return get_polymer_mass(mol, Descriptors.ExactMolWt)
    else:
        return Descriptors.ExactMolWt(mol)


def get_inchi_and_key(molfile):
    inchi = Chem.MolBlockToInchi(molfile)
    inchi_key = Chem.InchiToInchiKey(inchi)
    return inchi, inchi_key


def get_smiles(molfile):
    # we do apply RDKit chemistry model for SMILES generation
    mol = Chem.MolFromMolBlock(molfile)
    smiles = Chem.MolToSmiles(mol)
    return smiles


def get_polymer_formula(mol):
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

        formula = get_frag_formula(sub_mol)
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
    rest_formula = get_frag_formula(rwmol)

    if rest_formula:
        formulas.append(rest_formula)
    return ".".join(formulas)


def get_polymer_mass(mol, func):
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
    if rest_mass > 0.0: # remaining '*' have mass 0.0
        masses.append(str(rest_mass))
    return "+".join(masses)
