"""
Single Atom Analysis Module for MDemon

This module provides comprehensive single atom analysis capabilities including:
- Radial distribution function (RDF) analysis
- Diffusion coefficient calculation
- Coordination number analysis
- Parallel processing with automatic memory management

Classes:
    SingleAtomAnalyzer: Abstract base class for all single atom analyzers
    AnalysisConfig: Configuration management for analysis parameters
    AnalysisResult: Base class for analysis results

Functions:
    create_atom_selection_mask: Create atom selection masks with various criteria
    validate_universe: Validate Universe objects for single atom analysis
"""

from .base import (
    AnalysisConfig,
    AnalysisResult,
    SingleAtomAnalyzer,
    create_atom_selection_mask,
    validate_universe,
)
from .main import SingleAtomAnalysis
from .rdf import RDFAnalyzer, RDFResult

__all__ = [
    "SingleAtomAnalyzer",
    "AnalysisConfig",
    "AnalysisResult",
    "create_atom_selection_mask",
    "validate_universe",
    "RDFAnalyzer",
    "RDFResult",
    "SingleAtomAnalysis",
]

__version__ = "0.1.0"
__author__ = "MDemon Development Team"
