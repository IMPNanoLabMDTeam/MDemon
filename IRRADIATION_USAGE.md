# Waligorski-Zhang Calculator Usage Guide

The `WaligorskiZhangCalculator` has been integrated into the MDemon project and now uses the standardized constants system.

## Quick Start

### Basic Usage with Predefined Materials

```python
from MDemon.utils import WaligorskiZhangCalculator
from MDemon.constants import MATERIAL
import numpy as np

# Create calculator for β-Ga₂O₃
calc = WaligorskiZhangCalculator.from_material_preset('Ga2O3')

# Set up calculation parameters
radius_array = np.logspace(0, 3, 100)  # 1 to 1000 nm
ion_Z = 8  # Oxygen ion
ion_energy = 10.0  # MeV/amu
density = MATERIAL.get_density('Ga2O3_beta') * 1e-3  # Convert to g/cm³

# Calculate dose distribution
dose_array = np.array([
    calc.calculate_radial_dose(r, ion_Z, ion_energy, density)
    for r in radius_array
])

# Calculate temperature distribution
temperature_array = np.array([
    calc.dose_to_temperature(dose, density)
    for dose in dose_array
])
```

### Available Predefined Materials

```python
# List all available materials
materials = WaligorskiZhangCalculator.list_available_materials()
print(materials)
# Output: ['Ga2O3', 'Al2O3', 'SiO2', 'Si', 'GaN', 'water', 'SiC', 'diamond', 'polyethylene']
```

### Custom Material Properties

```python
# Define custom material (example: TiO2)
calc = WaligorskiZhangCalculator(
    molecular_weight=0.07987,  # kg/mol
    atoms_per_molecule=3,      # Ti + 2*O
    ionization_energy_eV=82.0  # eV
)

# Use the calculator
dose = calc.calculate_radial_dose(radius=10.0, Z=26, energy_MeV_per_amu=5.0, density_g_per_cm3=4.23)
temperature = calc.dose_to_temperature(dose, 4.23)
```

### Complete Analysis

```python
# Calculate complete radial energy distribution
results = calc.calculate_radial_energy_distribution(
    radius_array_nm=np.logspace(0, 2.5, 50),
    Z=40,  # Zr ion
    energy_MeV_per_amu=20.0,
    density_g_cm3=3.21
)

# Results contain:
# - radius: array of radial distances [nm]
# - dose_density: array of dose values [keV/nm³]
# - energy_per_molecule_eV: energy per molecule [eV]
# - energy_per_atom_eV: energy per atom [eV]
# - cumulative_energy: cumulative energy deposition
# - total_energy_keV: total deposited energy [keV]
```

## Key Changes from Original Code

1. **No more external constants dictionary**: All physical constants are now sourced from `MDemon.constants`
2. **No more target_material dictionary**: Material properties are passed directly to the constructor
3. **Predefined materials**: Common materials are available through `from_material_preset()`
4. **Integrated with MDemon**: The calculator is now part of the `MDemon.utils` module

## Constants Used

The calculator now uses the following constants from `MDemon.constants.IRRADIATION`:

- `electron_mass_keV`: ≈ 511.0 keV/c²
- `water_constant_keV_per_mm`: 1.15e-6 keV/mm
- `range_parameter_alpha`: 1.667 (for β > 0.03) or 1.079 (for β < 0.03)
- `eV_to_K_conversion`: ≈ 11604.5 K/eV
- `amu_to_MeV`: ≈ 931.5 MeV/c²
- `barkas_coefficient`: 125.0

## Example Script

See `examples/waligorski_zhang_example.py` for a complete working example with plotting capabilities.

## Material Properties Format

When defining custom materials, use the following units:
- `molecular_weight`: kg/mol
- `atoms_per_molecule`: dimensionless integer
- `ionization_energy_eV`: eV

## Integration with Other MDemon Features

The calculator can be used together with other MDemon modules:

```python
from MDemon.constants import MATERIAL, PHYSICS, CONVERSION
from MDemon.utils import WaligorskiZhangCalculator

# Get material density from MDemon database
density = MATERIAL.get_density('Ga2O3_beta') * 1e-3  # Convert kg/m³ to g/cm³

# Use conversion factors
energy_joules = 10.0 * CONVERSION.EV_TO_JOULE * 1e6  # 10 MeV to Joules

# Access fundamental constants
print(f"Avogadro constant: {PHYSICS.N_A:.2e} mol⁻¹")
``` 