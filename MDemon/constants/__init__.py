"""
Physical Constants for Molecular Dynamics Simulations

This module provides commonly used physical constants organized by category.
All constants are provided in SI units unless otherwise specified.

Usage:
    >>> from MDemon.constants import PHYSICS, CHEMISTRY, CONVERSION
    >>> from MDemon.constants import k_B, N_A, c_light
    >>> from MDemon.constants import ANGSTROM_TO_METER, EV_TO_JOULE

Categories:
    - PHYSICS: Fundamental physical constants
    - CHEMISTRY: Chemistry-related constants  
    - CONVERSION: Unit conversion factors
    - MATERIAL: Material properties constants
"""

from .fundamental import PHYSICS
from .chemistry import CHEMISTRY
from .conversion import CONVERSION
from .material import MATERIAL
from .irradiation import IRRADIATION
from . import utils

# 便捷导入：最常用的常量直接暴露在顶层
# 基本物理常量
k_B = PHYSICS.k_B  # Boltzmann constant
N_A = PHYSICS.N_A  # Avogadro's number
c_light = PHYSICS.c  # Speed of light
h_planck = PHYSICS.h  # Planck constant
e_charge = PHYSICS.e  # Elementary charge

# 常用转换因子
ANGSTROM_TO_METER = CONVERSION.ANGSTROM_TO_METER
EV_TO_JOULE = CONVERSION.EV_TO_JOULE
KCAL_TO_JOULE = CONVERSION.KCAL_TO_JOULE
HARTREE_TO_EV = CONVERSION.HARTREE_TO_EV

# 化学常量
R_gas = CHEMISTRY.R  # Gas constant
atomic_mass_unit = CHEMISTRY.u  # Atomic mass unit

__all__ = [
    'PHYSICS', 'CHEMISTRY', 'CONVERSION', 'MATERIAL', 'IRRADIATION',
    'k_B', 'N_A', 'c_light', 'h_planck', 'e_charge',
    'ANGSTROM_TO_METER', 'EV_TO_JOULE', 'KCAL_TO_JOULE', 'HARTREE_TO_EV',
    'R_gas', 'atomic_mass_unit'
] 