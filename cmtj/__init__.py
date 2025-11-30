"""
CMTJ - C Magnetic Tunnel Junctions

A Python package for magnetic tunnel junction simulations.
"""

from ._version import __version__

# Import all C++ extension functionality
try:
    import _cmtj
    from _cmtj import *
    # Also make sure we have the version available
    __all__ = ["__version__"] + [name for name in dir(_cmtj) if not name.startswith('_')]
except ImportError as e:
    import warnings
    warnings.warn(f"Could not import C++ extension: {e}")
    __all__ = ["__version__"]

# Import Python submodules - they will be available as cmtj.models, cmtj.utils, etc.
# The submodules are automatically available due to Python's package structure
