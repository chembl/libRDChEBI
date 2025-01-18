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
            self.run_test_for_data(data, get_monoisotopic_mass, "monoisotopic_mass", True)

    def test_avg_mass(self):
        for data in [r_group, m_r_groups, single_star, mixtures, atoms, isotopes]:
            self.run_test_for_data(data, get_avg_mass, "avg_mass", True)

    def test_mol_formula(self):
        for data in [r_group, m_r_groups, single_star, mixtures, atoms, isotopes]:
            self.run_test_for_data(data, get_small_molecule_formula, "mol_formula")

    def test_net_charge(self):
        for data in [r_group, m_r_groups, single_star, mixtures, atoms]:
            self.run_test_for_data(data, get_net_charge, "net_charge")


class TestPolymers(BaseChEBITest):
    """Tests for polymers and complex molecules"""

    def test_polymer_formula(self):
        for data in [polymers, extra_polymers]:
            self.run_test_for_data(data, get_polymer_formula, "mol_formula")

    def test_net_charge(self):
        for data in [polymers, extra_polymers]:
            self.run_test_for_data(data, get_net_charge, "net_charge")
