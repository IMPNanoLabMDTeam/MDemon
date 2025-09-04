"""
Single atom analysis for MDemon

This module provides various analysis capabilities for single atoms in molecular dynamics systems.

Main classes:
    SingleAtomAnalysis: Unified interface for all single atom analyses
    RDFAnalyzer: Radial distribution function analysis
    AngularAnalyzer: Angular distribution analysis
    CoordinationAnalyzer: Coordination number analysis
    RDFResult: Container for RDF analysis results
    AngularResult: Container for angular analysis results
    CoordinationResult: Container for coordination analysis results
"""

from .angular import AngularAnalyzer, AngularResult
from .base import AnalysisConfig, AnalysisResult, SingleAtomAnalyzer
from .coordination import CoordinationAnalyzer, CoordinationResult
from .main import SingleAtomAnalysis
from .rdf import RDFAnalyzer, RDFResult

__all__ = [
    # Main interface
    "SingleAtomAnalysis",
    # Base classes
    "SingleAtomAnalyzer",
    "AnalysisResult",
    "AnalysisConfig",
    # RDF analysis
    "RDFAnalyzer",
    "RDFResult",
    # Angular analysis
    "AngularAnalyzer",
    "AngularResult",
    # Coordination analysis
    "CoordinationAnalyzer",
    "CoordinationResult",
]

# Version information
__version__ = "1.0.0"

# Analysis types
AVAILABLE_ANALYSES = {
    "rdf": "Radial Distribution Function",
    "angular": "Angular Distribution Function",
    "coordination": "Coordination Number Analysis",
    "diffusion": "Diffusion Analysis (not implemented)",
}


def get_available_analyses():
    """
    Get information about available analysis types

    Returns
    -------
    dict
        Dictionary of available analysis types and descriptions
    """
    return AVAILABLE_ANALYSES.copy()
