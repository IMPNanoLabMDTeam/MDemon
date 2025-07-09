"""
Irradiation and Radiation Physics Constants

Constants specifically used for radiation physics calculations,
particularly for the Waligorski-Zhang radial dose distribution model.
"""

from .fundamental import PHYSICS
from .conversion import CONVERSION


class IrradiationConstants:
    """Container for irradiation physics constants"""
    
    # Electron mass in keV/c² (for relativistic calculations)
    electron_mass_keV = PHYSICS.m_e * PHYSICS.c**2 / CONVERSION.EV_TO_JOULE / 1000.0  # ≈ 511.0 keV/c²
    
    # Water constant for dose calculations [keV/mm]
    # This is an empirical constant from Waligorski-Zhang model
    water_constant_keV_per_mm = 8.5  # keV/mm
    
    # Range parameter alpha for Waligorski-Zhang model
    # Depends on target material and ion velocity
    range_parameter_alpha = 1.667  # Dimensionless (for β > 0.03)
    range_parameter_alpha_slow = 1.079  # Dimensionless (for β < 0.03)
    
    # ionization energy in Waligorski-Zhang model
    ionization_energy_eV = 10.0   # eV

    # WZ range parameter
    wz_range_g_per_cm2_keV_per_alpha = 6e-6  # g/cm²/keV**α

    # Beta velocity threshold for range parameter selection
    beta_threshold = 0.03  # Dimensionless
    
    # Barkas correction coefficient
    barkas_coefficient = 125.0  # For effective charge calculation
    
    # eV to Kelvin conversion factor (energy per particle to temperature)
    # Derived from k_B: E[eV] * eV_to_K = T[K]
    eV_to_K_conversion = CONVERSION.EV_TO_JOULE / PHYSICS.k_B  # ≈ 11604.5 K/eV
    
    # Atomic mass unit in MeV/c²
    amu_to_MeV = PHYSICS.u * PHYSICS.c**2 / (CONVERSION.EV_TO_JOULE * 1e6)  # ≈ 931.5 MeV/c²
    
    @classmethod
    def get_range_parameter_alpha(cls, beta):
        """Get appropriate range parameter alpha based on ion velocity"""
        return cls.range_parameter_alpha_slow if beta < cls.beta_threshold else cls.range_parameter_alpha
    
    @classmethod
    def electron_mass_MeV(cls):
        """Get electron mass in MeV/c²"""
        return cls.electron_mass_keV / 1000.0
    
    @classmethod
    def list_constants(cls):
        """Return a dictionary of all constants with descriptions"""
        descriptions = {
            'electron_mass_keV': 'Electron rest mass [keV/c²]',
            'water_constant_keV_per_mm': 'Water constant for Waligorski-Zhang model [keV/mm]',
            'range_parameter_alpha': 'Range parameter for fast ions (β > 0.03) [dimensionless]',
            'range_parameter_alpha_slow': 'Range parameter for slow ions (β < 0.03) [dimensionless]',
            'beta_threshold': 'Velocity threshold for range parameter selection [dimensionless]',
            'barkas_coefficient': 'Barkas correction coefficient [dimensionless]',
            'eV_to_K_conversion': 'Energy per particle to temperature conversion [K/eV]',
            'amu_to_MeV': 'Atomic mass unit in energy units [MeV/c²]'
        }
        
        constants = {}
        for name in descriptions:
            if hasattr(cls, name):
                constants[name] = {
                    'value': getattr(cls, name),
                    'description': descriptions[name]
                }
        return constants


# Create the singleton instance
IRRADIATION = IrradiationConstants() 