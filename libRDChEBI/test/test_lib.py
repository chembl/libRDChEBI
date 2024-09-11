from ..descriptors import (
    get_net_charge,
    get_small_molecule_formula,
    get_polymer_formula,
    get_avg_mass,
    get_monoisotopic_mass,
)
from .mols import (
    r_group,
    m_r_groups,
    single_star,
    polymers,
    mixtures,
    atoms,
    extra_polymers,
    isotopes,
)
from pytest import approx


class TestRGroupMols:
    def test_monoisotopic_mass(self):
        for key, mol in r_group.items():
            assert get_monoisotopic_mass(mol["molfile"]) == approx(
                mol["monoisotopic_mass"]
            ), f"ChEBI:{key}"

    def test_avgMass(self):
        for key, mol in r_group.items():
            assert get_avg_mass(mol["molfile"]) == approx(
                mol["avg_mass"]
            ), f"ChEBI:{key}"

    def test_molFormula(self):
        for key, mol in r_group.items():
            assert (
                get_small_molecule_formula(mol["molfile"]) == mol["mol_formula"]
            ), f"ChEBI:{key}"

    def test_netCharge(self):
        for key, mol in r_group.items():
            assert get_net_charge(mol["molfile"]) == mol["net_charge"], f"ChEBI:{key}"


class TestMultipleRGroupMols:
    def test_monoisotopic_mass(self):
        for key, mol in m_r_groups.items():
            assert get_monoisotopic_mass(mol["molfile"]) == approx(
                mol["monoisotopic_mass"]
            ), f"ChEBI:{key}"

    def test_avgMass(self):
        for key, mol in m_r_groups.items():
            assert get_avg_mass(mol["molfile"]) == approx(
                mol["avg_mass"]
            ), f"ChEBI:{key}"

    def test_molFormula(self):
        for key, mol in m_r_groups.items():
            assert (
                get_small_molecule_formula(mol["molfile"]) == mol["mol_formula"]
            ), f"ChEBI:{key}"

    def test_netCharge(self):
        for key, mol in m_r_groups.items():
            assert get_net_charge(mol["molfile"]) == mol["net_charge"], f"ChEBI:{key}"


class TestSingleStar:
    def test_monoisotopic_mass(self):
        for key, mol in single_star.items():
            assert get_monoisotopic_mass(mol["molfile"]) == approx(
                mol["monoisotopic_mass"]
            ), f"ChEBI:{key}"

    def test_avgMass(self):
        for key, mol in single_star.items():
            assert get_avg_mass(mol["molfile"]) == approx(
                mol["avg_mass"]
            ), f"ChEBI:{key}"

    def test_molFormula(self):
        for key, mol in single_star.items():
            assert (
                get_small_molecule_formula(mol["molfile"]) == mol["mol_formula"]
            ), f"ChEBI:{key}"

    def test_netCharge(self):
        for key, mol in single_star.items():
            assert get_net_charge(mol["molfile"]) == mol["net_charge"], f"ChEBI:{key}"


# class TestPolymers:

#     def test_molFormula(self):
#         for key, mol in polymers.items():
#             if mol["mol_formula"] is not None:
#                 assert (
#                     get_polymer_formula(mol["molfile"]) == mol["mol_formula"]
#                 ), f"ChEBI:{key}"

#     def test_netCharge(self):
#         for key, mol in polymers.items():
#             if mol["net_charge"] is not None:
#                 assert (
#                     get_net_charge(mol["molfile"]) == mol["net_charge"]
#                 ), f"ChEBI:{key}"


class TestMixtures:
    def test_monoisotopic_mass(self):
        for key, mol in mixtures.items():
            if mol["monoisotopic_mass"] is not None:
                assert get_monoisotopic_mass(mol["molfile"]) == approx(
                    mol["monoisotopic_mass"]
                ), f"ChEBI:{key}"

    def test_avgMass(self):
        for key, mol in mixtures.items():
            if mol["avg_mass"] is not None:
                assert get_avg_mass(mol["molfile"]) == approx(
                    mol["avg_mass"]
                ), f"ChEBI:{key}"

    def test_molFormula(self):
        for key, mol in mixtures.items():
            if mol["mol_formula"] is not None:
                assert (
                    get_small_molecule_formula(mol["molfile"]) == mol["mol_formula"]
                ), f"ChEBI:{key}"

    def test_netCharge(self):
        for key, mol in mixtures.items():
            if mol["net_charge"] is not None:
                assert (
                    get_net_charge(mol["molfile"]) == mol["net_charge"]
                ), f"ChEBI:{key}"


class TestAtoms:
    def test_monoisotopic_mass(self):
        for key, mol in atoms.items():
            if mol["monoisotopic_mass"] is not None:
                assert get_monoisotopic_mass(mol["molfile"]) == approx(
                    mol["monoisotopic_mass"]
                ), f"ChEBI:{key}"

    def test_avgMass(self):
        for key, mol in atoms.items():
            if mol["avg_mass"] is not None:
                assert get_avg_mass(mol["molfile"]) == approx(
                    mol["avg_mass"]
                ), f"ChEBI:{key}"

    def test_molFormula(self):
        for key, mol in atoms.items():
            if mol["mol_formula"] is not None:
                assert (
                    get_small_molecule_formula(mol["molfile"]) == mol["mol_formula"]
                ), f"ChEBI:{key}"

    def test_netCharge(self):
        for key, mol in atoms.items():
            if mol["net_charge"] is not None:
                assert (
                    get_net_charge(mol["molfile"]) == mol["net_charge"]
                ), f"ChEBI:{key}"


# class TestExtraPolymers:

#     def test_molFormula(self):
#         for key, mol in extra_polymers.items():
#             if mol["mol_formula"] is not None:
#                 assert (
#                     get_polymer_formula(mol["molfile"]) == mol["mol_formula"]
#                 ), f"ChEBI:{key}"

#     def test_netCharge(self):
#         for key, mol in extra_polymers.items():
#             if mol["net_charge"] is not None:
#                 assert (
#                     get_net_charge(mol["molfile"]) == mol["net_charge"]
#                 ), f"ChEBI:{key}"


class TestIsotopes:
    def test_monoisotopic_mass(self):
        for key, mol in isotopes.items():
            if mol["monoisotopic_mass"] is not None:
                assert get_monoisotopic_mass(mol["molfile"]) == approx(
                    mol["monoisotopic_mass"]
                ), f"ChEBI:{key}"

    def test_avgMass(self):
        for key, mol in isotopes.items():
            if mol["avg_mass"] is not None:
                assert get_avg_mass(mol["molfile"]) == approx(
                    mol["avg_mass"]
                ), f"ChEBI:{key}"

    def test_molFormula(self):
        for key, mol in isotopes.items():
            if mol["mol_formula"] is not None:
                assert (
                    get_small_molecule_formula(mol["molfile"]) == mol["mol_formula"]
                ), f"ChEBI:{key}"
