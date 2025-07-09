"""
MDemon Physical Constants Usage Examples

This script demonstrates how to use the physical constants module in MDemon
for common molecular dynamics simulation calculations.
"""

# Import MDemon constants
from MDemon.constants import PHYSICS, CHEMISTRY, CONVERSION, MATERIAL
from MDemon.constants import k_B, N_A, c_light, ANGSTROM_TO_METER, EV_TO_JOULE
from MDemon.constants import utils


def demonstrate_basic_constants():
    """Demonstrate basic physical constants usage"""
    print("=== Basic Physical Constants ===")
    print(f"Boltzmann constant: {k_B:.6e} J/K")
    print(f"Avogadro's number: {N_A:.6e} mol⁻¹")
    print(f"Speed of light: {c_light:.0f} m/s")
    print(f"Elementary charge: {PHYSICS.e:.6e} C")
    print(f"Planck constant: {PHYSICS.h:.6e} J⋅s")
    print()


def demonstrate_unit_conversions():
    """Demonstrate unit conversion capabilities"""
    print("=== Unit Conversions ===")
    
    # Energy conversions
    energy_ev = 1.0  # 1 eV
    energy_j = utils.convert_energy(energy_ev, 'eV', 'J')
    energy_kcal = utils.convert_energy(energy_ev, 'eV', 'kcal/mol')
    print(f"1 eV = {energy_j:.6e} J = {energy_kcal:.4f} kcal/mol")
    
    # Length conversions  
    length_angstrom = 5.0  # 5 Å
    length_nm = utils.convert_length(length_angstrom, 'Å', 'nm')
    length_m = utils.convert_length(length_angstrom, 'Å', 'm')
    print(f"5 Å = {length_nm:.2f} nm = {length_m:.2e} m")
    
    # Temperature conversions
    temp_c = 25.0  # 25°C
    temp_k = utils.convert_temperature(temp_c, 'C', 'K')
    temp_f = utils.convert_temperature(temp_c, 'C', 'F')
    print(f"25°C = {temp_k:.2f} K = {temp_f:.1f}°F")
    print()


def demonstrate_thermal_calculations():
    """Demonstrate thermal energy calculations"""
    print("=== Thermal Energy Calculations ===")
    
    # Room temperature
    temp = 298.15  # K
    
    # Thermal energy k_B * T
    thermal_energy_j = utils.boltzmann_energy(temp, 'J')
    thermal_energy_ev = utils.boltzmann_energy(temp, 'eV')
    thermal_energy_kcal = utils.boltzmann_energy(temp, 'kcal/mol')
    
    print(f"Thermal energy at {temp} K:")
    print(f"  {thermal_energy_j:.6e} J")
    print(f"  {thermal_energy_ev:.4f} eV")
    print(f"  {thermal_energy_kcal:.4f} kcal/mol")
    
    # Thermal velocity for different atoms
    masses = {
        'H': 1.008 * CHEMISTRY.u,  # Hydrogen
        'C': 12.011 * CHEMISTRY.u,  # Carbon
        'O': 15.999 * CHEMISTRY.u,  # Oxygen
        'Ga': 69.723 * CHEMISTRY.u,  # Gallium
    }
    
    print(f"\nThermal velocities at {temp} K:")
    for atom, mass in masses.items():
        v_thermal = utils.thermal_velocity(mass, temp)
        print(f"  {atom}: {v_thermal:.0f} m/s")
    print()


def demonstrate_material_properties():
    """Demonstrate material property lookups"""
    print("=== Material Properties ===")
    
    # Properties for Ga₂O₃
    materials = ['Ga2O3_beta', 'Al2O3', 'Si', 'Cu']
    properties = ['density', 'melting_point', 'thermal_conductivity']
    
    for material in materials:
        print(f"\n{material}:")
        for prop in properties:
            value = utils.get_material_property(material, prop)
            if value is not None:
                units = {
                    'density': 'kg/m³',
                    'melting_point': 'K',
                    'thermal_conductivity': 'W/(m⋅K)'
                }
                print(f"  {prop}: {value} {units[prop]}")
    
    # Lattice parameters
    print(f"\nLattice parameters:")
    ga2o3_params = MATERIAL.get_lattice_parameters('Ga2O3_beta')
    if ga2o3_params:
        print(f"β-Ga₂O₃: a={ga2o3_params['a']:.3f} Å, b={ga2o3_params['b']:.3f} Å, "
              f"c={ga2o3_params['c']:.3f} Å, β={ga2o3_params['beta']:.2f}°")
    print()


def demonstrate_chemistry_constants():
    """Demonstrate chemistry-related constants"""
    print("=== Chemistry Constants ===")
    
    # Gas constant in different units
    print(f"Gas constant R:")
    print(f"  {CHEMISTRY.R:.4f} J/(mol⋅K)")
    print(f"  {CHEMISTRY.R_cal:.4f} cal/(mol⋅K)")
    print(f"  {CHEMISTRY.R_liter_atm:.5f} L⋅atm/(mol⋅K)")
    
    # Standard conditions
    print(f"\nStandard conditions:")
    print(f"  STP temperature: {CHEMISTRY.STP_temperature} K")
    print(f"  STP pressure: {CHEMISTRY.STP_pressure} Pa")
    print(f"  Molar volume at STP: {CHEMISTRY.STP_molar_volume:.4f} m³/mol")
    
    # Bond lengths and van der Waals radii
    print(f"\nTypical bond lengths:")
    bonds = ['C-C', 'C=C', 'C-H', 'O-H']
    for bond in bonds:
        length = CHEMISTRY.get_bond_length(bond)
        if length:
            length_angstrom = utils.convert_length(length, 'm', 'Å')
            print(f"  {bond}: {length_angstrom:.2f} Å")
    
    print(f"\nVan der Waals radii:")
    elements = ['H', 'C', 'N', 'O']
    for element in elements:
        radius = CHEMISTRY.get_vdw_radius(element)
        if radius:
            radius_angstrom = utils.convert_length(radius, 'm', 'Å')
            print(f"  {element}: {radius_angstrom:.2f} Å")
    print()


def demonstrate_advanced_calculations():
    """Demonstrate more advanced physics calculations"""
    print("=== Advanced Calculations ===")
    
    # Debye screening length in electrolyte
    temp = 298.15  # K
    ionic_strength = 0.1  # mol/L (typical for physiological conditions)
    debye_len_m = utils.debye_length(temp, ionic_strength, 'm')
    debye_len_nm = utils.debye_length(temp, ionic_strength, 'nm')
    print(f"Debye length at {temp} K, I={ionic_strength} M:")
    print(f"  {debye_len_m:.2e} m = {debye_len_nm:.1f} nm")
    
    # Frequency from energy
    photon_energy_ev = 2.5  # eV (visible light)
    photon_energy_j = utils.convert_energy(photon_energy_ev, 'eV', 'J')
    frequency_hz = utils.planck_frequency(photon_energy_j, 'Hz')
    frequency_thz = utils.planck_frequency(photon_energy_j, 'THz')
    print(f"\nPhoton with energy {photon_energy_ev} eV:")
    print(f"  Frequency: {frequency_hz:.2e} Hz = {frequency_thz:.1f} THz")
    
    # Thermal de Broglie wavelength
    mass_electron = PHYSICS.m_e
    temp = 298.15  # K
    import math
    h = PHYSICS.h
    k_B = PHYSICS.k_B
    thermal_wavelength = h / math.sqrt(2 * math.pi * mass_electron * k_B * temp)
    thermal_wavelength_angstrom = utils.convert_length(thermal_wavelength, 'm', 'Å')
    print(f"\nThermal de Broglie wavelength of electron at {temp} K:")
    print(f"  {thermal_wavelength_angstrom:.2f} Å")
    print()


def main():
    """Run all demonstration functions"""
    print("MDemon Physical Constants Usage Examples\n")
    print("=" * 50)
    
    demonstrate_basic_constants()
    demonstrate_unit_conversions()
    demonstrate_thermal_calculations()
    demonstrate_material_properties()
    demonstrate_chemistry_constants()
    demonstrate_advanced_calculations()
    
    # Show available units
    print("=== Available Units ===")
    available_units = utils.list_available_units()
    for category, units in available_units.items():
        print(f"{category.capitalize()}: {', '.join(units)}")
    
    print("\n" + "=" * 50)
    print("All examples completed successfully!")


if __name__ == "__main__":
    main() 