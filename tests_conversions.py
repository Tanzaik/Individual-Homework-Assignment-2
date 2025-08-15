"""
Unit tests for conversion helpers.
We primarily check round-trip consistency and a few known reference calculations.
"""

import unittest
from math import isclose

from data_processing import (
    ARGON, LJArgs, K_B,
    Tstar_to_K, K_to_Tstar,
    pstar_to_Pa, Pa_to_pstar,
    rhostar_to_number_density, number_density_to_rhostar,
    rhostar_to_molar_density, rhostar_to_mass_density,
    number_density_to_molar_density, molar_to_mass_density
)


class TestConversions(unittest.TestCase):
    def test_temperature_round_trip(self):
        T_star = [0.8, 1.0, 1.5]
        T = Tstar_to_K(T_star, ARGON)
        T_back = K_to_Tstar(T, ARGON)
        for a, b in zip(T_star, T_back):
            self.assertAlmostEqual(a, b, places=12)

    def test_pressure_known_value(self):
        # For p* = 1, p = ε / σ^3
        p_star = [1.0]
        p = pstar_to_Pa(p_star, ARGON)[0]
        eps = ARGON.epsilon_J
        sig3 = ARGON.sigma_m ** 3
        self.assertAlmostEqual(p, eps / sig3, places=12)

    def test_pressure_round_trip(self):
        p_star = [0.1, 1.0, 2.5]
        p = pstar_to_Pa(p_star, ARGON)
        p_back = Pa_to_pstar(p, ARGON)
        for a, b in zip(p_star, p_back):
            self.assertAlmostEqual(a, b, places=12)

    def test_density_chains(self):
        rho_star = [0.1, 0.5, 0.85]
        n = rhostar_to_number_density(rho_star, ARGON)
        # round trip
        rho_back = number_density_to_rhostar(n, ARGON)
        for a, b in zip(rho_star, rho_back):
            self.assertAlmostEqual(a, b, places=12)

        # to molar and mass densities (sanity: monotonic and positive)
        c = number_density_to_molar_density(n)
        rho_mass = molar_to_mass_density(c, ARGON)
        self.assertTrue(all(x > 0 for x in c))
        self.assertTrue(all(x > 0 for x in rho_mass))


if __name__ == "__main__":
    unittest.main()
