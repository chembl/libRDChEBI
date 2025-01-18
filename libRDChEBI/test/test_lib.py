from ..descriptors import (
    get_net_charge,
    get_molformula,
    get_avg_mass,
    get_monoisotopic_mass,
    get_mass_from_formula,
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


class BaseChEBITest:
    """Base class for ChEBI test cases"""

    def run_test_for_data(self, data, func, property_name, approx_test=False):
        """Generic test runner for molecule data"""
        for key, mol in data.items():
            expected = mol[property_name]
            if expected is not None:
                if approx_test:
                    assert func(mol["molfile"]) == approx(expected), f"ChEBI:{key}"
                else:
                    assert func(mol["molfile"]) == expected, f"ChEBI:{key}"


class TestSmallMolecules(BaseChEBITest):
    """Tests for regular molecules"""

    def test_monoisotopic_mass(self):
        for data in [r_group, m_r_groups, single_star, mixtures, atoms, isotopes]:
            self.run_test_for_data(
                data, get_monoisotopic_mass, "monoisotopic_mass", True
            )

    def test_avg_mass(self):
        for data in [r_group, m_r_groups, single_star, mixtures, atoms, isotopes]:
            self.run_test_for_data(data, get_avg_mass, "avg_mass", True)

    def test_mol_formula(self):
        for data in [r_group, m_r_groups, single_star, mixtures, atoms, isotopes]:
            self.run_test_for_data(data, get_molformula, "mol_formula")

    def test_net_charge(self):
        for data in [r_group, m_r_groups, single_star, mixtures, atoms]:
            self.run_test_for_data(data, get_net_charge, "net_charge")


class TestPolymers(BaseChEBITest):
    """Tests for polymers and complex molecules"""

    def test_polymer_formula(self):
        for data in [polymers, extra_polymers]:
            self.run_test_for_data(data, get_molformula, "mol_formula")

    def test_net_charge(self):
        for data in [polymers, extra_polymers]:
            self.run_test_for_data(data, get_net_charge, "net_charge")


class TestMassFromFormula:
    """Tests for calculating mass from molecular formula"""

    def test_mass_from_formula(self):
        # average mass (average=True)
        assert get_mass_from_formula("H2O", average=True) == approx(18.015)
        assert get_mass_from_formula("CH4", average=True) == approx(16.043)
        assert get_mass_from_formula("NaCl", average=True) == approx(58.443)
        assert get_mass_from_formula("C6H6", average=True) == approx(78.114)
        assert get_mass_from_formula("R2", average=True) == approx(0.0)  # R groups have no mass

        # monoisotopic mass (average=False)
        assert get_mass_from_formula("H2O", average=False) == approx(18.010565)
        assert get_mass_from_formula("CH4", average=False) == approx(16.031300)
        assert get_mass_from_formula("NaCl", average=False) == approx(57.958621)
        assert get_mass_from_formula("C6H6", average=False) == approx(78.046950)
        assert get_mass_from_formula("R2", average=False) == approx(0.0)  # R groups have no mass

        # invalid formulas
        assert get_mass_from_formula("XxYy", average=True) is None
        assert get_mass_from_formula("ABC", average=False) is None
