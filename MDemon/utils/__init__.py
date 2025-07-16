"""
Utilities for Molecular Dynamics Analysis

This module contains utility classes and functions for various MD analysis tasks,
including radiation physics calculations, Dask-based parallel processing, and
spatial subdivision for efficient neighbor searching.

Available utilities:
    - WaligorskiZhangCalculator: Radial dose distribution calculator for ion irradiation
    - DaskParallelManager: Dask-based parallel processing with automatic memory management
    - SpatialSubdivision: Spatial subdivision for efficient neighbor searching
"""

from .irradiation import WaligorskiZhangCalculator
from .parallel import DaskParallelManager
from .spatial import SpatialSubdivision

__all__ = [
    "WaligorskiZhangCalculator",
    "DaskParallelManager",
    "SpatialSubdivision",
]
