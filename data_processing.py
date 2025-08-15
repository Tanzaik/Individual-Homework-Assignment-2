"""
Data parsing + reduced-unit <-> real-unit conversions for a Lennard–Jones (LJ) fluid.
Target species: Argon (default constants provided, but can be overridden).
"""

from __future__ import annotations
import csv
from dataclasses import dataclass
from typing import List, Tuple, Iterable, Optional

# ---- Physical constants ----
K_B = 1.380649e-23        # Boltzmann constant [J/K]
N_A = 6.02214076e23       # Avogadro's number [1/mol]

@dataclass(frozen=True)
class LJArgs:
    """Parameters needed to convert LJ reduced units <-> real units for a monatomic gas."""
    sigma_angstrom: float         # LJ size parameter σ in Å
    epsilon_over_k_B_K: float     # ε/k_B in K
    molar_mass_kg_per_mol: float  # molar mass (Argon ~ 0.039948 kg/mol)

    @property
    def sigma_m(self) -> float:
        return self.sigma_angstrom * 1e-10

    @property
    def epsilon_J(self) -> float:
        # ε = k_B * (ε/k_B)
        return K_B * self.epsilon_over_k_B_K


# Default: Argon
ARGON = LJArgs(sigma_angstrom=3.4, epsilon_over_k_B_K=120.0, molar_mass_kg_per_mol=0.039948)


# ---------------- CSV parsing ----------------
def parse_monte_carlo_data(filepath: str) -> Tuple[List[float], List[float], List[float]]:
    """
    Parse reduced LJ data from a CSV: [T*, rho*, p*, ...]
    Expected columns (at least): 
        col0=T* (reduced temperature), col1=rho* (reduced density), col2.. maybe others, col4 could be p*
    Because formats vary, we try p* in a few common columns.
    Returns: (T_star, rho_star, p_star)
    """
    T_star, rho_star, p_star = [], [], []
    with open(filepath, "r", newline="") as f:
        reader = csv.reader(f)
        header = next(reader, None)  # skip header if present
        for row in reader:
            if not row or all(x.strip()=="" for x in row):
                continue
            T_star.append(float(row[0]))
            rho_star.append(float(row[1]))
            # Try common positions for reduced pressure
            p_val = None
            for idx in (2, 3, 4, 5):
                if idx < len(row):
                    try:
                        p_val = float(row[idx])
                        break
                    except ValueError:
                        continue
            if p_val is None:
                raise ValueError("Could not locate reduced pressure column (p*) in Monte Carlo CSV.")
            p_star.append(p_val)
    return T_star, rho_star, p_star


def parse_experimental_data(filepath: str) -> Tuple[List[float], List[float], List[float]]:
    """
    Parse experimental Argon data CSV.
    Expected columns: T [K], density [kg/m^3 or mol/m^3], pressure [Pa or bar]
    We only read the first three numeric columns and leave units handling to the caller.
    Returns: (T, density, pressure)
    """
    T, rho, p = [], [], []
    with open(filepath, "r", newline="") as f:
        reader = csv.reader(f)
        header = next(reader, None)  # skip header
        for row in reader:
            if not row or all(x.strip()=="" for x in row):
                continue
            T.append(float(row[0]))
            rho.append(float(row[1]))
            p.append(float(row[2]))
    return T, rho, p


# ---------------- Conversions (reduced <-> real) ----------------
# Definitions used (standard for LJ):
#  - Reduced temperature:   T* = (k_B * T) / ε
#  - Reduced pressure:      p* = p * σ^3 / ε
#  - Reduced number density rho* = n * σ^3, where n is number density [1/m^3]
#
# From these:
#  - T  [K]      = T* * ε / k_B = T* * (ε/k_B)
#  - p  [Pa]     = p* * ε / σ^3
#  - n  [1/m^3]  = rho* / σ^3
#  - molar density c [mol/m^3] = n / N_A
#  - mass density ρ_mass [kg/m^3] = c * molar_mass


def Tstar_to_K(T_star: Iterable[float], params: LJArgs = ARGON) -> List[float]:
    return [tstar * params.epsilon_over_k_B_K for tstar in T_star]


def K_to_Tstar(T_K: Iterable[float], params: LJArgs = ARGON) -> List[float]:
    return [t / params.epsilon_over_k_B_K for t in T_K]


def pstar_to_Pa(p_star: Iterable[float], params: LJArgs = ARGON) -> List[float]:
    sig3 = params.sigma_m ** 3
    eps = params.epsilon_J
    return [ps * eps / sig3 for ps in p_star]


def Pa_to_pstar(p_Pa: Iterable[float], params: LJArgs = ARGON) -> List[float]:
    sig3 = params.sigma_m ** 3
    eps = params.epsilon_J
    return [p * sig3 / eps for p in p_Pa]


def rhostar_to_number_density(rho_star: Iterable[float], params: LJArgs = ARGON) -> List[float]:
    sig3 = params.sigma_m ** 3
    return [rs / sig3 for rs in rho_star]  # [1/m^3]


def number_density_to_rhostar(n_per_m3: Iterable[float], params: LJArgs = ARGON) -> List[float]:
    sig3 = params.sigma_m ** 3
    return [n * sig3 for n in n_per_m3]


def number_density_to_molar_density(n_per_m3: Iterable[float]) -> List[float]:
    return [n / N_A for n in n_per_m3]  # [mol/m^3]


def molar_density_to_number_density(c_mol_per_m3: Iterable[float]) -> List[float]:
    return [c * N_A for c in c_mol_per_m3]  # [1/m^3]


def molar_to_mass_density(c_mol_per_m3: Iterable[float], params: LJArgs = ARGON) -> List[float]:
    M = params.molar_mass_kg_per_mol
    return [c * M for c in c_mol_per_m3]  # [kg/m^3]


def rhostar_to_mass_density(rho_star: Iterable[float], params: LJArgs = ARGON) -> List[float]:
    n = rhostar_to_number_density(rho_star, params)
    c = number_density_to_molar_density(n)
    return molar_to_mass_density(c, params)


def rhostar_to_molar_density(rho_star: Iterable[float], params: LJArgs = ARGON) -> List[float]:
    n = rhostar_to_number_density(rho_star, params)
    return number_density_to_molar_density(n)


# ---------------- Convenience batch conversion ----------------
def convert_reduced_to_real(
    T_star: Iterable[float],
    rho_star: Iterable[float],
    p_star: Iterable[float],
    params: LJArgs = ARGON,
    pressure_unit: str = "Pa",
    density_unit: str = "kg/m^3",
) -> Tuple[List[float], List[float], List[float]]:
    """
    Convert arrays of (T*, rho*, p*) to (T [K], rho [density_unit], p [pressure_unit]).
    pressure_unit in {"Pa", "bar", "MPa"}, density_unit in {"kg/m^3", "mol/m^3", "1/m^3"}
    """
    T_K = Tstar_to_K(T_star, params)
    p_Pa = pstar_to_Pa(p_star, params)

    # Pressure unit handling
    if pressure_unit == "Pa":
        p_out = p_Pa
    elif pressure_unit == "bar":
        p_out = [p / 1e5 for p in p_Pa]
    elif pressure_unit == "MPa":
        p_out = [p / 1e6 for p in p_Pa]
    else:
        raise ValueError("Unsupported pressure_unit")

    # Density unit handling
    if density_unit == "kg/m^3":
        rho_out = rhostar_to_mass_density(rho_star, params)
    elif density_unit == "mol/m^3":
        rho_out = rhostar_to_molar_density(rho_star, params)
    elif density_unit == "1/m^3":
        rho_out = rhostar_to_number_density(rho_star, params)
    else:
        raise ValueError("Unsupported density_unit")

    return T_K, rho_out, p_out
