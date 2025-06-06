"""
Enhancement Engine - Comprehensive genetic enhancement simulation and analysis.

This package provides tools for:
- CRISPR guide design and optimization
- Off-target prediction and safety assessment  
- Enhancement effect simulation
- Comprehensive reporting and analysis

Author: Tim Röhr
Version: 0.1.0
"""

__version__ = "0.1.0"
__author__ = "Tim Röhr"
__email__ = "tim.roehr@outlook.com"

# Core imports for easy access
from .core.engine import EnhancementEngine
from .models.data_classes import (
    GeneInfo,
    GuideRNA, 
    EnhancementReport,
    SafetyScore,
    VariantEffect
)

# Version info
VERSION_INFO = {
    "major": 0,
    "minor": 1,
    "patch": 0,
    "release": "alpha"
}

# Default configuration
DEFAULT_CONFIG = {
    "ncbi_email": None,  # Must be set by user
    "cache_enabled": True,
    "cache_dir": "data/cached_sequences",
    "log_level": "INFO",
    "default_organism": "homo sapiens",
    "default_cas_type": "cas9"
}

__all__ = [
    "EnhancementEngine",
    "GeneInfo", 
    "GuideRNA",
    "EnhancementReport",
    "SafetyScore",
    "VariantEffect",
    "__version__",
    "VERSION_INFO",
    "DEFAULT_CONFIG"
]