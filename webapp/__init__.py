"""
Flask web interface for the Enhancement Engine.

This module provides both a simple form-based interface and a JSON API
for genetic enhancement analysis.

Features:
- Web form for gene analysis
- JSON API endpoints
- Error handling and logging
- Support for both local development and container deployment
"""

import os
import sys
from pathlib import Path

# Add parent directory to path for local development
current_dir = Path(__file__).parent
parent_dir = current_dir.parent
if str(parent_dir) not in sys.path:
    sys.path.insert(0, str(parent_dir))

from flask import Flask

# Version info
__version__ = "0.1.0"
__author__ = "Tim Röhr"

# Import main app factory
try:
    from .run import create_app, main
    __all__ = ['create_app', 'main', '__version__', '__author__']
except ImportError:
    # Fallback for development
    __all__ = ['__version__', '__author__']

# Default configuration
DEFAULT_CONFIG = {
    'HOST': '127.0.0.1',
    'PORT': 5000,
    'DEBUG': True,
    'SECRET_KEY': 'development-key',
    'DEMO_EMAIL': 'demo@example.com'
}