"""
Utilities for Molecular Dynamics Analysis

This module contains utility classes and functions for various MD analysis tasks,
including radiation physics calculations and other specialized computations.

Available utilities:
    - WaligorskiZhangCalculator: Radial dose distribution calculator for ion irradiation
    - MaterialPresets: Predefined material properties for common target materials
"""

from .irradiation import WaligorskiZhangCalculator

__all__ = [
    'WaligorskiZhangCalculator',
] 