"""
Unit Conversion Factors

Conversion factors for common units used in molecular dynamics simulations.
All factors are multiplicative: quantity_in_SI = quantity_in_unit * factor
"""

import math


class ConversionFactors:
    """Container for unit conversion factors"""
    
    # Length conversions (to meters)
    ANGSTROM_TO_METER = 1e-10
    NANOMETER_TO_METER = 1e-9
    PICOMETER_TO_METER = 1e-12
    MICROMETER_TO_METER = 1e-6
    CENTIMETER_TO_METER = 1e-2
    MILLIMETER_TO_METER = 1e-3
    BOHR_TO_METER = 5.29177210903e-11
    
    # Inverse length conversions (to m⁻¹)
    METER_TO_ANGSTROM = 1e10
    METER_TO_NANOMETER = 1e9
    METER_TO_PICOMETER = 1e12
    METER_TO_BOHR = 1.8897259886e10
    
    # Energy conversions (to Joules)
    EV_TO_JOULE = 1.602176634e-19
    KCAL_TO_JOULE = 4184.0
    CAL_TO_JOULE = 4.184
    HARTREE_TO_JOULE = 4.3597447222071e-18
    RY_TO_JOULE = 2.1798723611035e-18  # Rydberg
    WAVENUMBER_TO_JOULE = 1.98630e-23  # cm⁻¹ to J
    
    # Energy conversions (commonly used in computational chemistry)
    HARTREE_TO_EV = 27.211386245988
    HARTREE_TO_KCAL = 627.5094740631
    EV_TO_KCAL = 23.060541
    EV_TO_WAVENUMBER = 8065.5439  # eV to cm⁻¹
    KCAL_TO_EV = 0.04336414
    
    # Temperature conversions
    CELSIUS_TO_KELVIN = 273.15  # Add this to Celsius to get Kelvin
    FAHRENHEIT_TO_KELVIN_FACTOR = 5.0/9.0  # Multiply (F-32) by this
    RANKINE_TO_KELVIN = 5.0/9.0  # Multiply Rankine by this
    
    # Pressure conversions (to Pascal)
    ATM_TO_PASCAL = 101325.0
    BAR_TO_PASCAL = 100000.0
    TORR_TO_PASCAL = 133.322368
    PSI_TO_PASCAL = 6894.757
    MMHG_TO_PASCAL = 133.322368
    
    # Time conversions (to seconds)
    FEMTOSECOND_TO_SECOND = 1e-15
    PICOSECOND_TO_SECOND = 1e-12
    NANOSECOND_TO_SECOND = 1e-9
    MICROSECOND_TO_SECOND = 1e-6
    MILLISECOND_TO_SECOND = 1e-3
    MINUTE_TO_SECOND = 60.0
    HOUR_TO_SECOND = 3600.0
    DAY_TO_SECOND = 86400.0
    
    # Mass conversions (to kg)
    GRAM_TO_KG = 1e-3
    AMU_TO_KG = 1.66053906660e-27
    POUND_TO_KG = 0.45359237
    
    # Force conversions (to Newton)
    DYNE_TO_NEWTON = 1e-5
    POUND_FORCE_TO_NEWTON = 4.448222
    
    # Volume conversions (to m³)
    LITER_TO_CUBIC_METER = 1e-3
    CUBIC_CM_TO_CUBIC_METER = 1e-6
    CUBIC_ANGSTROM_TO_CUBIC_METER = 1e-30
    
    # Angle conversions (to radians)
    DEGREE_TO_RADIAN = math.pi / 180.0
    RADIAN_TO_DEGREE = 180.0 / math.pi
    
    # Electric field conversions (to V/m)
    ATOMIC_UNIT_FIELD_TO_VM = 5.1422067476e11  # E_h/(e⋅a₀) to V/m
    
    # Dipole moment conversions (to C⋅m)
    DEBYE_TO_CM = 3.33564095e-30
    ATOMIC_UNIT_DIPOLE_TO_CM = 8.4783536255e-30  # e⋅a₀ to C⋅m
    
    # Frequency conversions (to Hz)
    WAVENUMBER_TO_HZ = 2.99792458e10  # cm⁻¹ to Hz
    THZ_TO_HZ = 1e12
    GHZ_TO_HZ = 1e9
    MHZ_TO_HZ = 1e6
    
    @classmethod
    def temperature_celsius_to_kelvin(cls, celsius):
        """Convert Celsius to Kelvin"""
        return celsius + cls.CELSIUS_TO_KELVIN
    
    @classmethod
    def temperature_fahrenheit_to_kelvin(cls, fahrenheit):
        """Convert Fahrenheit to Kelvin"""
        return (fahrenheit - 32.0) * cls.FAHRENHEIT_TO_KELVIN_FACTOR + cls.CELSIUS_TO_KELVIN
    
    @classmethod
    def energy_units_available(cls):
        """Return available energy conversion factors"""
        return {
            'eV': cls.EV_TO_JOULE,
            'kcal/mol': cls.KCAL_TO_JOULE,
            'cal': cls.CAL_TO_JOULE,
            'Hartree': cls.HARTREE_TO_JOULE,
            'Rydberg': cls.RY_TO_JOULE,
            'cm⁻¹': cls.WAVENUMBER_TO_JOULE
        }
    
    @classmethod
    def length_units_available(cls):
        """Return available length conversion factors"""
        return {
            'Å': cls.ANGSTROM_TO_METER,
            'nm': cls.NANOMETER_TO_METER,
            'pm': cls.PICOMETER_TO_METER,
            'μm': cls.MICROMETER_TO_METER,
            'cm': cls.CENTIMETER_TO_METER,
            'mm': cls.MILLIMETER_TO_METER,
            'bohr': cls.BOHR_TO_METER
        }


# Create the singleton instance
CONVERSION = ConversionFactors() 