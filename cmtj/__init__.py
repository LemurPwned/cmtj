"""
CMTJ - C Magnetic Tunnel Junctions

A Python package for magnetic tunnel junction simulations.
"""

from ._version import __version__

# Import the C++ extension module packaged inside ``cmtj``
try:
    from . import _cmtj
except ImportError as e:
    import warnings
    warnings.warn(f"Could not import C++ extension: {e}")
else:
    for name in dir(_cmtj):
        if not name.startswith("_"):
            globals()[name] = getattr(_cmtj, name)

__all__ = ["__version__"]
try:
    __all__.extend([n for n in dir(_cmtj) if not n.startswith("_")])
except Exception:
    pass
