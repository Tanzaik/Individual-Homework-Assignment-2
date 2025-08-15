# Individual-Homework-Assignment-2

This program converts the Monte Carlo simulation results given in *reduced* Lennard–Jones (LJ) units to **real** SI units for **Argon**, and compare to experimental measurements at a common temperature by plotting

## The files in this repository consist of:

- `data_processing.py` — parsing functions for the provided CSV files and **conversion helpers** that map reduced LJ quantities to real units for Argon.
- `main.py` — a runnable script that reads the reduced data, converts to real units, optionally overlays experimental data, and generates a plot.
- `tests_conversions.py` — **unit tests** that cover each conversion function.
- (You provide) `nist_data.csv` — Monte Carlo reduced data with at least `T*`, `rho*`, and `p*` columns.
- (Optional) `nist_tabulated_argon.csv` — real Argon data (e.g., at 108 K) for comparison.


## Background — reduced vs. real units

For a Lennard–Jones fluid, reduced (dimensionless) variables are defined by the species-specific parameters `σ` (size) and `ε` (well depth). Using the standard LJ definitions:

- Reduced temperature:  \u200b`T* = (k_B * T) / ε`  
- Reduced number density: `ρ* = n * σ^3` where `n` is number density `[1/m^3]`  
- Reduced pressure: `p* = p * σ^3 / ε`

From these we convert to real units (for Argon, defaults used: σ = 3.4 Å, ε/k_B = 120 K, M = 39.948 g/mol):

- `T [K] = T* * (ε/k_B)`  
- `p [Pa] = p* * ε / σ^3`  
- `n [1/m^3] = ρ* / σ^3`  
- `c [mol/m^3] = n / N_A`  
- `ρ_mass [kg/m^3] = c * M`

These are implemented in `data_processing.py` and unit-tested.


## How to run

1. **Install requirements** (only matplotlib for plotting; csv is stdlib):
   ```bash
   pip install matplotlib
   ```

2. **Place the input files** at the project root:
   - `nist_data.csv` — reduced LJ data (first two columns must be `T*`, `ρ*`; pressure can be any of the next columns).
   - Optionally `nist_tabulated_argon.csv` — experimental Argon data with three columns `[T, density, pressure]` in whatever units that file uses.

3. **Convert + Plot** (example):
   ```bash
   python main.py \
     --monte-carlo-csv nist_data.csv \
     --exp-csv nist_tabulated_argon.csv \
     --density-unit "kg/m^3" \
     --pressure-unit bar \
     --out-csv converted_argon.csv \
     --plot-path argon_plot.png
   ```

   - `--density-unit` options: `kg/m^3` (mass density), `mol/m^3`, or raw number density `1/m^3`  
   - `--pressure-unit` options: `Pa`, `bar`, `MPa`

4. **Run tests**:
   ```bash
   python -m unittest tests_conversions.py -v
   ```



## Reflection

During the conversion, the most error-prone steps are unit discipline and interpreting what experimental “density” represents. LJ reduced density corresponds to **number density**, so to compare with typical experimental tables (often given as **mass density** in `kg/m^3`), we map  
`ρ* → n [1/m^3] → c [mol/m^3] → ρ_mass [kg/m^3]` using Avogadro’s number and Argon’s molar mass.  
Pressure conversion follows directly from `p = p* ε / σ^3`. Guarding against common pitfalls, the code centralizes constants (σ, ε/k_B, M) in a dataclass, so swapping species is trivial and tests verify round-trip consistency.


## File format assumptions

- **Monte Carlo CSV** must have `T*` in column 0 and `ρ*` in column 1; `p*` can be in column 2–5 (script tries these automatically).  
- **Experimental CSV** is read as `[T, density, pressure]` without unit conversion because sources vary; pass flags that match its units when plotting for apples-to-apples overlays.
