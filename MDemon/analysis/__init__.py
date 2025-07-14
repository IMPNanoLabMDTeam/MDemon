"""
Analysis Module for MDemon

This module provides various analysis capabilities for molecular dynamics data,
including single atom analysis, structural analysis, and thermodynamic analysis.

Submodules:
    single_atom: Single atom analysis including RDF, diffusion, and coordination analysis
"""

from . import single_atom

__all__ = [
    "single_atom",
]

__version__ = "0.1.0"
__author__ = "MDemon Development Team"
