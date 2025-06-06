"""
Data models and structures for Enhancement Engine.

This module contains all data classes, constants, and type definitions
used throughout the Enhancement Engine package.
"""

from .data_classes import *
from .constants import *

__all__ = [
    # Data classes
    "GeneInfo",
    "GuideRNA",
    "EnhancementReport",
    "SafetyScore",
    "VariantEffect",
    "OffTarget",
    "PAMSite",
    "EfficiencyScore",
    "EnhancementGain",
    "VariantPosition",
    "CodingRegion",
    # Constants
    "ENHANCEMENT_GENES",
    "PAM_PATTERNS",
    "CAS_TYPES",
    "NCBI_DATABASES",
    "DEFAULT_PARAMETERS",
]
