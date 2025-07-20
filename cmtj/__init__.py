"""
CMTJ - C Magnetic Tunnel Junctions

A Python package for magnetic tunnel junction simulations.
"""

from ._version import __version__

# Import the C++ extension module
try:
    # Import the compiled C++ extension (renamed to _cmtj to avoid conflicts)
    from .. import _cmtj
    
    # Import all public symbols from the C++ extension
    for name in dir(_cmtj):
        if not name.startswith('_'):
            globals()[name] = getattr(_cmtj, name)
            
except ImportError:
    try:
        # Fallback: try direct import for development/testing
        import _cmtj
        
        # Import all public symbols from the C++ extension
        for name in dir(_cmtj):
            if not name.startswith('_'):
                globals()[name] = getattr(_cmtj, name)
                
    except ImportError as e:
        import warnings
        warnings.warn(f"Could not import C++ extension: {e}")

__all__ = ["__version__"]
