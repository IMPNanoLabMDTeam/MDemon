"""
Chemistry-related Physical Constants

Constants commonly used in molecular dynamics simulations and computational chemistry.
Units are SI unless otherwise specified.
"""

from .fundamental import PHYSICS


class ChemistryConstants:
    """Container for chemistry-related constants"""
    
    # Gas constants
    R = PHYSICS.R  # Molar gas constant [J/(mol⋅K)]
    R_cal = 1.987204144  # Gas constant [cal/(mol⋅K)]
    R_liter_atm = 0.08205736  # Gas constant [L⋅atm/(mol⋅K)]
    
    # Atomic and molecular masses
    u = PHYSICS.u  # Atomic mass unit [kg]
    u_g = 1.66053906660e-24  # Atomic mass unit [g]
    
    # Standard conditions
    STP_temperature = 273.15  # Standard temperature [K]
    STP_pressure = 101325.0  # Standard pressure [Pa]
    STP_molar_volume = 0.022413969545052  # Molar volume at STP [m³/mol]
    
    # Faraday constant
    F = PHYSICS.e * PHYSICS.N_A  # Faraday constant [C/mol]
    
    # Standard state conditions
    standard_pressure = 100000.0  # Standard pressure (IUPAC) [Pa]
    standard_temperature = 298.15  # Standard temperature [K]
    
    # Calorie conversion
    cal_to_joule = 4.184  # Thermochemical calorie [J/cal]
    
    # Common molecular properties
    water_density_stp = 999.84  # Water density at STP [kg/m³]
    water_molar_mass = 0.01801528  # Water molar mass [kg/mol]
    
    # Bond lengths (typical values in meters)
    bond_lengths = {
        'C-C': 1.54e-10,  # Single bond [m]
        'C=C': 1.34e-10,  # Double bond [m] 
        'C≡C': 1.20e-10,  # Triple bond [m]
        'C-H': 1.09e-10,  # C-H bond [m]
        'O-H': 0.96e-10,  # O-H bond [m]
        'N-H': 1.01e-10,  # N-H bond [m]
        'C-O': 1.43e-10,  # C-O single bond [m]
        'C=O': 1.23e-10,  # C=O double bond [m]
        'C-N': 1.47e-10,  # C-N single bond [m]
        'C=N': 1.29e-10,  # C=N double bond [m]
        'N-N': 1.45e-10,  # N-N single bond [m]
        'N=N': 1.25e-10,  # N=N double bond [m]
        'O-O': 1.48e-10,  # O-O single bond [m]
    }
    
    # Van der Waals radii (in meters)
    vdw_radii = {
        'H': 1.20e-10,
        'C': 1.70e-10,
        'N': 1.55e-10,
        'O': 1.52e-10,
        'F': 1.47e-10,
        'P': 1.80e-10,
        'S': 1.80e-10,
        'Cl': 1.75e-10,
        'Ar': 1.88e-10,
        'Br': 1.85e-10,
        'I': 1.98e-10,
    }
    
    # Electronegativity (Pauling scale)
    electronegativity = {
        'H': 2.20,
        'C': 2.55,
        'N': 3.04,
        'O': 3.44,
        'F': 3.98,
        'Na': 0.93,
        'Mg': 1.31,
        'Al': 1.61,
        'Si': 1.90,
        'P': 2.19,
        'S': 2.58,
        'Cl': 3.16,
        'K': 0.82,
        'Ca': 1.00,
        'Br': 2.96,
        'I': 2.66,
    }
    
    @classmethod
    def get_bond_length(cls, bond_type):
        """Get typical bond length for a given bond type"""
        return cls.bond_lengths.get(bond_type, None)
    
    @classmethod
    def get_vdw_radius(cls, element):
        """Get van der Waals radius for an element"""
        return cls.vdw_radii.get(element, None)
    
    @classmethod
    def get_electronegativity(cls, element):
        """Get electronegativity for an element"""
        return cls.electronegativity.get(element, None)


# Create the singleton instance
CHEMISTRY = ChemistryConstants() 