"""
Example script demonstrating energy conservation verification
in the Waligorski-Zhang calculator.

This script shows how the radial integral of the normalized dose
should equal the input energy loss multiplied by the g-factor.
"""

import numpy as np

from MDemon.utils.irradiation import WaligorskiZhangCalculator

# Example 1: Xenon ion in PET (Polyethylene Terephthalate)
print("=" * 70)
print("Example 1: Xe ion in PET")
print("=" * 70)

calculator_pet = WaligorskiZhangCalculator(
    molecular_weight=192.17,  # PET molecular weight (g/mol)
    atoms_per_molecule=22,  # C10H8O4 has 22 atoms
    density_g_per_cm3=1.38,  # PET density
    ion_name="Xe",
    ion_energy_MeV_per_amu=1.0,  # 1 MeV/amu
    ion_Z=54,  # Xenon atomic number
    energy_loss_keV_per_um=10.5,  # Energy loss in PET
    g_factor=0.17,
)

# Create radial grid
radius_array = np.linspace(0.1, 50, 500)  # 0.1 to 50 nm with 500 points

# Calculate energy distribution
results_pet = calculator_pet.calculate_radial_energy_distribution(radius_array)

print(
    f"\nRadial grid: {len(radius_array)} points from {radius_array[0]:.1f} to {radius_array[-1]:.1f} nm"
)
print(f"Grid spacing: average Δr = {np.mean(np.diff(radius_array)):.3f} nm")

# Example 2: Gallium ion in water
print("\n" + "=" * 70)
print("Example 2: Ga ion in H2O (for comparison)")
print("=" * 70)

calculator_water = WaligorskiZhangCalculator(
    molecular_weight=18.015,  # Water molecular weight
    atoms_per_molecule=3,  # H2O
    density_g_per_cm3=1.0,  # Water density
    ion_name="Ga",
    ion_energy_MeV_per_amu=2.0,
    ion_Z=31,  # Gallium atomic number
    energy_loss_keV_per_um=8.2,
    g_factor=0.2,
)

# Create a finer radial grid for better accuracy
radius_array_fine = np.linspace(0.05, 100, 1000)  # Finer grid

results_water = calculator_water.calculate_radial_energy_distribution(radius_array_fine)

print(
    f"\nRadial grid: {len(radius_array_fine)} points from {radius_array_fine[0]:.2f} to {radius_array_fine[-1]:.1f} nm"
)
print(f"Grid spacing: average Δr = {np.mean(np.diff(radius_array_fine)):.3f} nm")

# Example 3: Test grid refinement effect
print("\n" + "=" * 70)
print("Example 3: Grid Refinement Test (same system, different grids)")
print("=" * 70)

test_grids = [
    ("Coarse", np.linspace(0.1, 50, 100)),
    ("Medium", np.linspace(0.1, 50, 500)),
    ("Fine", np.linspace(0.1, 50, 2000)),
]

calculator_test = WaligorskiZhangCalculator(
    molecular_weight=192.17,
    atoms_per_molecule=22,
    density_g_per_cm3=1.38,
    ion_name="Xe",
    ion_energy_MeV_per_amu=1.0,
    ion_Z=54,
    energy_loss_keV_per_um=10.5,
    g_factor=0.17,
)

print("\nGrid refinement comparison:")
print("-" * 70)
for grid_name, grid in test_grids:
    results = calculator_test.calculate_radial_energy_distribution(grid)
    print(
        f"{grid_name:10s} grid ({len(grid):4d} points): Error = {results['energy_conservation_error_percent']:.4f}%"
    )

print("\n" + "=" * 70)
print("Summary:")
print("=" * 70)
print("The energy conservation verification ensures that:")
print("  ∫₀^∞ D(r) · 2πr dr = (dE/dx) × g")
print("where:")
print("  D(r) = normalized dose density [keV/nm³]")
print("  r = radial distance [nm]")
print("  dE/dx = linear energy transfer [keV/nm]")
print("  g = g-factor (energy deposition efficiency)")
print("\nFiner radial grids generally provide better energy conservation.")
print("=" * 70)
