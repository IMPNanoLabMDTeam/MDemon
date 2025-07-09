"""
Utility functions for working with physical constants and unit conversions.

This module provides convenient functions for common operations in molecular dynamics simulations.
"""

from .fundamental import PHYSICS
from .chemistry import CHEMISTRY
from .conversion import CONVERSION
from .material import MATERIAL


def convert_energy(value, from_unit, to_unit='J'):
    """
    Convert energy between different units.
    
    Parameters:
    -----------
    value : float
        Energy value to convert
    from_unit : str
        Source unit ('J', 'eV', 'kcal/mol', 'cal', 'Hartree', 'Ry', 'cm-1')
    to_unit : str
        Target unit (default: 'J')
    
    Returns:
    --------
    float
        Converted energy value
    """
    # Unit mapping
    unit_map = {
        'J': 1.0,
        'eV': CONVERSION.EV_TO_JOULE,
        'kcal/mol': CONVERSION.KCAL_TO_JOULE,
        'kcal': CONVERSION.KCAL_TO_JOULE,
        'cal': CONVERSION.CAL_TO_JOULE,
        'Hartree': CONVERSION.HARTREE_TO_JOULE,
        'Ha': CONVERSION.HARTREE_TO_JOULE,
        'Ry': CONVERSION.RY_TO_JOULE,
        'cm-1': CONVERSION.WAVENUMBER_TO_JOULE,
        'wavenumber': CONVERSION.WAVENUMBER_TO_JOULE,
    }
    
    if from_unit not in unit_map:
        raise ValueError(f"Unknown source unit: {from_unit}")
    if to_unit not in unit_map:
        raise ValueError(f"Unknown target unit: {to_unit}")
    
    # Convert to Joules first, then to target unit
    value_in_joules = value * unit_map[from_unit]
    converted_value = value_in_joules / unit_map[to_unit]
    
    return converted_value


def convert_length(value, from_unit, to_unit='m'):
    """
    Convert length between different units.
    
    Parameters:
    -----------
    value : float
        Length value to convert
    from_unit : str
        Source unit ('m', 'Å', 'nm', 'pm', 'μm', 'cm', 'mm', 'bohr')
    to_unit : str
        Target unit (default: 'm')
    
    Returns:
    --------
    float
        Converted length value
    """
    unit_map = {
        'm': 1.0,
        'Å': CONVERSION.ANGSTROM_TO_METER,
        'angstrom': CONVERSION.ANGSTROM_TO_METER,
        'nm': CONVERSION.NANOMETER_TO_METER,
        'pm': CONVERSION.PICOMETER_TO_METER,
        'μm': CONVERSION.MICROMETER_TO_METER,
        'um': CONVERSION.MICROMETER_TO_METER,
        'cm': CONVERSION.CENTIMETER_TO_METER,
        'mm': CONVERSION.MILLIMETER_TO_METER,
        'bohr': CONVERSION.BOHR_TO_METER,
    }
    
    if from_unit not in unit_map:
        raise ValueError(f"Unknown source unit: {from_unit}")
    if to_unit not in unit_map:
        raise ValueError(f"Unknown target unit: {to_unit}")
    
    # Convert to meters first, then to target unit
    value_in_meters = value * unit_map[from_unit]
    converted_value = value_in_meters / unit_map[to_unit]
    
    return converted_value


def convert_temperature(value, from_unit, to_unit='K'):
    """
    Convert temperature between different units.
    
    Parameters:
    -----------
    value : float
        Temperature value to convert
    from_unit : str
        Source unit ('K', 'C', 'F', 'R')
    to_unit : str
        Target unit (default: 'K')
    
    Returns:
    --------
    float
        Converted temperature value
    """
    # First convert to Kelvin
    if from_unit == 'K':
        kelvin_value = value
    elif from_unit == 'C':
        kelvin_value = CONVERSION.temperature_celsius_to_kelvin(value)
    elif from_unit == 'F':
        kelvin_value = CONVERSION.temperature_fahrenheit_to_kelvin(value)
    elif from_unit == 'R':  # Rankine
        kelvin_value = value * CONVERSION.RANKINE_TO_KELVIN
    else:
        raise ValueError(f"Unknown source unit: {from_unit}")
    
    # Then convert from Kelvin to target unit
    if to_unit == 'K':
        return kelvin_value
    elif to_unit == 'C':
        return kelvin_value - CONVERSION.CELSIUS_TO_KELVIN
    elif to_unit == 'F':
        return (kelvin_value - CONVERSION.CELSIUS_TO_KELVIN) / CONVERSION.FAHRENHEIT_TO_KELVIN_FACTOR + 32.0
    elif to_unit == 'R':
        return kelvin_value / CONVERSION.RANKINE_TO_KELVIN
    else:
        raise ValueError(f"Unknown target unit: {to_unit}")


def boltzmann_energy(temperature, unit='J'):
    """
    Calculate thermal energy k_B * T.
    
    Parameters:
    -----------
    temperature : float
        Temperature in Kelvin
    unit : str
        Energy unit for result (default: 'J')
    
    Returns:
    --------
    float
        Thermal energy k_B * T
    """
    thermal_energy_j = PHYSICS.k_B * temperature
    return convert_energy(thermal_energy_j, 'J', unit)


def debye_length(temperature, ionic_strength, unit='m'):
    """
    Calculate Debye screening length.
    
    Parameters:
    -----------
    temperature : float
        Temperature in Kelvin
    ionic_strength : float
        Ionic strength in mol/L
    unit : str
        Length unit for result (default: 'm')
    
    Returns:
    --------
    float
        Debye length
    """
    import math
    
    epsilon_r = 78.4  # Relative permittivity of water at 25°C
    epsilon = epsilon_r * PHYSICS.epsilon_0
    k_b_T = PHYSICS.k_B * temperature
    e = PHYSICS.e
    N_A = PHYSICS.N_A
    
    # Convert ionic strength from mol/L to m⁻³
    ionic_strength_m3 = ionic_strength * 1000 * N_A
    
    # Debye length formula
    debye_length_m = math.sqrt(epsilon * k_b_T / (2 * e**2 * ionic_strength_m3))
    
    return convert_length(debye_length_m, 'm', unit)


def thermal_velocity(mass, temperature, unit='m/s'):
    """
    Calculate thermal velocity (root mean square velocity).
    
    Parameters:
    -----------
    mass : float
        Particle mass in kg
    temperature : float
        Temperature in Kelvin
    unit : str
        Velocity unit for result (default: 'm/s')
    
    Returns:
    --------
    float
        RMS thermal velocity
    """
    import math
    
    k_b_T = PHYSICS.k_B * temperature
    v_thermal = math.sqrt(3 * k_b_T / mass)
    
    # For now, just return in m/s (could extend for other velocity units)
    if unit != 'm/s':
        raise NotImplementedError("Only m/s velocity unit supported currently")
    
    return v_thermal


def planck_frequency(energy, unit='Hz'):
    """
    Calculate frequency from energy using Planck relation E = h*ν.
    
    Parameters:
    -----------
    energy : float
        Energy in Joules
    unit : str
        Frequency unit for result (default: 'Hz')
    
    Returns:
    --------
    float
        Frequency
    """
    frequency_hz = energy / PHYSICS.h
    
    if unit == 'Hz':
        return frequency_hz
    elif unit == 'THz':
        return frequency_hz / CONVERSION.THZ_TO_HZ
    elif unit == 'GHz':
        return frequency_hz / CONVERSION.GHZ_TO_HZ
    elif unit == 'MHz':
        return frequency_hz / CONVERSION.MHZ_TO_HZ
    else:
        raise ValueError(f"Unknown frequency unit: {unit}")


def get_material_property(material, property_name):
    """
    Get a material property by name.
    
    Parameters:
    -----------
    material : str
        Material name
    property_name : str
        Property name ('density', 'melting_point', 'thermal_conductivity', 
                      'bulk_modulus', 'youngs_modulus')
    
    Returns:
    --------
    float or None
        Property value or None if not found
    """
    property_methods = {
        'density': MATERIAL.get_density,
        'melting_point': MATERIAL.get_melting_point,
        'thermal_conductivity': MATERIAL.get_thermal_conductivity,
        'bulk_modulus': MATERIAL.get_bulk_modulus,
        'youngs_modulus': MATERIAL.get_youngs_modulus,
    }
    
    if property_name not in property_methods:
        available = list(property_methods.keys())
        raise ValueError(f"Unknown property: {property_name}. Available: {available}")
    
    return property_methods[property_name](material)


def list_available_units():
    """
    List all available units for conversions.
    
    Returns:
    --------
    dict
        Dictionary with unit categories and available units
    """
    return {
        'energy': ['J', 'eV', 'kcal/mol', 'cal', 'Hartree', 'Ry', 'cm-1'],
        'length': ['m', 'Å', 'nm', 'pm', 'μm', 'cm', 'mm', 'bohr'],
        'temperature': ['K', 'C', 'F', 'R'],
        'frequency': ['Hz', 'THz', 'GHz', 'MHz'],
    } 