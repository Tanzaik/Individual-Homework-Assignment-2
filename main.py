# Individual Homework Assignment 2
# Name: Tanzim Amin  
# Date: August 15, 2025

"""
Run end-to-end:
1) Parse reduced NIST Monte Carlo data
2) Convert to real units for Argon
3) (Optional) Parse experimental data for comparison
4) Plot Pressure vs Density at a given isotherm
"""

from pathlib import Path
import argparse
import csv

import matplotlib.pyplot as plt

from data_processing import (
    parse_monte_carlo_data, parse_experimental_data,
    convert_reduced_to_real, ARGON
)


def write_converted_csv(outpath, T, rho, p, density_unit, pressure_unit):
    outpath = Path(outpath)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    with open(outpath, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([f"T [K]", f"density [{density_unit}]", f"pressure [{pressure_unit}]"])
        for ti, ri, pi in zip(T, rho, p):
            w.writerow([ti, ri, pi])


def main(args):
    # 1) Parse reduced data
    T_star, rho_star, p_star = parse_monte_carlo_data(args.monte_carlo_csv)

    # 2) Convert to real units
    T, rho, p = convert_reduced_to_real(
        T_star, rho_star, p_star, ARGON,
        pressure_unit=args.pressure_unit,
        density_unit=args.density_unit
    )

    # 3) Optionally save converted table
    if args.out_csv:
        write_converted_csv(args.out_csv, T, rho, p, args.density_unit, args.pressure_unit)

    # 4) Optional experimental comparison: if provided, we just plot on same axes
    plt.figure()
    plt.xlabel(f"Density [{args.density_unit}]")
    plt.ylabel(f"Pressure [{args.pressure_unit}]")
    plt.title("Argon LJ Monte Carlo vs (optional) Experiment")

    plt.scatter(rho, p, label="Monte Carlo (converted)")

    if args.exp_csv:
        exp_T, exp_rho, exp_p = parse_experimental_data(args.exp_csv)
        # No unit conversions attempted here because source units vary by dataset.
        # If your file is in bar and kg/m^3, pass --pressure-unit bar --density-unit "kg/m^3"
        plt.scatter(exp_rho, exp_p, marker="x", label="Experiment (as-is)")

    plt.legend()
    plt.tight_layout()

    if args.plot_path:
        Path(args.plot_path).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(args.plot_path, dpi=200)
    else:
        plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert LJ reduced data to real units for Argon and visualize.")
    parser.add_argument("--monte-carlo-csv", required=True, help="Path to reduced-unit CSV, e.g., nist_data.csv")
    parser.add_argument("--exp-csv", default=None, help="Optional experimental CSV for overlay")
    parser.add_argument("--density-unit", default="kg/m^3", choices=["kg/m^3", "mol/m^3", "1/m^3"], help="Output density unit")
    parser.add_argument("--pressure-unit", default="bar", choices=["Pa", "bar", "MPa"], help="Output pressure unit")
    parser.add_argument("--out-csv", default=None, help="Optional path to save converted table")
    parser.add_argument("--plot-path", default=None, help="Optional path to save plot (PNG)")
    main(parser.parse_args())
