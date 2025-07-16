"""
Selection utilities for MDemon

This package provides various selection tools for filtering atoms based on
geometric and spatial criteria.

Classes:
    CylindricalSelector: Select atoms within a cylindrical region
"""

from .cylindrical import CylindricalSelector

__all__ = ["CylindricalSelector"]
