"""
Utilities for Molecular Dynamics Analysis

This module contains utility classes and functions for various MD analysis tasks,
including radiation physics calculations and Dask-based parallel processing.

Available utilities:
    - WaligorskiZhangCalculator: Radial dose distribution calculator for ion irradiation
    - DaskParallelManager: Dask-based parallel processing with automatic memory management
"""

from .irradiation import WaligorskiZhangCalculator
from .parallel import DaskParallelManager

__all__ = [
    "WaligorskiZhangCalculator",
    "DaskParallelManager",
]
