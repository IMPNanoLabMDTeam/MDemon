"""
Fundamental Physical Constants

All values are based on 2018 CODATA internationally recommended values.
Units are SI unless otherwise specified.

Reference: https://physics.nist.gov/cuu/Constants/
"""

import math


class PhysicsConstants:
    """Container for fundamental physical constants"""
    
    # Universal constants
    c = 299792458.0  # Speed of light in vacuum [m/s]
    h = 6.62607015e-34  # Planck constant [J⋅s]
    hbar = h / (2 * math.pi)  # Reduced Planck constant [J⋅s]
    
    # Electromagnetic constants  
    e = 1.602176634e-19  # Elementary charge [C]
    epsilon_0 = 8.8541878128e-12  # Vacuum permittivity [F/m]
    mu_0 = 1.25663706212e-6  # Vacuum permeability [H/m]
    
    # Thermodynamic constants
    k_B = 1.380649e-23  # Boltzmann constant [J/K]
    N_A = 6.02214076e23  # Avogadro constant [mol⁻¹]
    R = k_B * N_A  # Molar gas constant [J/(mol⋅K)]
    
    # Particle masses
    m_e = 9.1093837015e-31  # Electron mass [kg]
    m_p = 1.67262192369e-27  # Proton mass [kg]
    m_n = 1.67492749804e-27  # Neutron mass [kg]
    u = 1.66053906660e-27  # Atomic mass unit [kg]
    
    # Gravitational constant
    G = 6.67430e-11  # Gravitational constant [m³/(kg⋅s²)]
    
    # Fine structure constant
    alpha = 7.2973525693e-3  # Fine-structure constant [dimensionless]
    
    # Rydberg constant
    R_inf = 10973731.568160  # Rydberg constant [m⁻¹]
    
    # Bohr radius
    a_0 = 5.29177210903e-11  # Bohr radius [m]
    
    # Classical electron radius
    r_e = 2.8179403262e-15  # Classical electron radius [m]
    
    @classmethod
    def list_constants(cls):
        """Return a dictionary of all constants with descriptions"""
        descriptions = {
            'c': 'Speed of light in vacuum [m/s]',
            'h': 'Planck constant [J⋅s]',
            'hbar': 'Reduced Planck constant [J⋅s]',
            'e': 'Elementary charge [C]',
            'epsilon_0': 'Vacuum permittivity [F/m]',
            'mu_0': 'Vacuum permeability [H/m]',
            'k_B': 'Boltzmann constant [J/K]',
            'N_A': 'Avogadro constant [mol⁻¹]',
            'R': 'Molar gas constant [J/(mol⋅K)]',
            'm_e': 'Electron mass [kg]',
            'm_p': 'Proton mass [kg]',
            'm_n': 'Neutron mass [kg]',
            'u': 'Atomic mass unit [kg]',
            'G': 'Gravitational constant [m³/(kg⋅s²)]',
            'alpha': 'Fine-structure constant [dimensionless]',
            'R_inf': 'Rydberg constant [m⁻¹]',
            'a_0': 'Bohr radius [m]',
            'r_e': 'Classical electron radius [m]'
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
PHYSICS = PhysicsConstants() 