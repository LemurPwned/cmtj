"""
Unified constants module for CMTJ.

This module provides a Python interface to the C++ constants system,
allowing runtime modification of physical constants used throughout
the CMTJ library.

Usage:
    from cmtj.utils.constants import Constants

    # Get current values
    mu0 = Constants.mu0()
    gamma = Constants.gamma()

    # Modify constants at runtime
    Constants.set_mu0(12.57e-7)
    Constants.set_gamma(28024e6)

    # Reset to defaults
    Constants.reset_to_defaults()
"""

import math
import warnings

# Try to import the C++ constants module, fall back to pure Python if not available
try:
    from .. import constants as _cpp_constants

    _HAS_CPP_CONSTANTS = True
except ImportError:
    _HAS_CPP_CONSTANTS = False
    warnings.warn(
        "C++ constants module not available. Using pure Python constants.",
        ImportWarning,
        stacklevel=2,
    )


class Constants:
    """
    Unified constants interface that delegates to C++ when available,
    falls back to Python constants otherwise.
    """

    # Fallback Python constants (used when C++ module is not available)
    _FALLBACK_CONSTANTS = {
        "mu0": 12.566e-7,  # magnetic permeability of free space
        "gyromagnetic_ratio": 2.211e5,  # m/As
        "gamma": 28024e6,  # Hz/T
        "gamma_rad": 1.76e11,  # rad / (s * T)
        "TtoAm": 795774.715459,
        "AmtoT": 1.0 / 795774.715459,
        "OetoAm": 79.57747,
        "AmtoOe": 1.0 / 79.57747,
        "hplanck": 6.6260e-34,
        "hbar": 6.6260e-34 / (2 * math.pi),
        "echarge": -1.602e-19,
        "me": 9.109e-31,
        "boltzmann": 1.380649e-23,
    }

    # Storage for pure Python mode
    _python_constants = _FALLBACK_CONSTANTS.copy()

    @classmethod
    def _use_cpp(cls) -> bool:
        """Check if C++ constants should be used."""
        return _HAS_CPP_CONSTANTS

    # Magnetic permeability
    @classmethod
    def mu0(cls) -> float:
        """Get magnetic permeability of free space [H/m]."""
        if cls._use_cpp():
            return _cpp_constants.PhysicalConstants.magnetic_permeability()
        return cls._python_constants["mu0"]

    @classmethod
    def set_mu0(cls, value: float) -> None:
        """Set magnetic permeability of free space [H/m]."""
        if cls._use_cpp():
            _cpp_constants.PhysicalConstants.set_magnetic_permeability(value)
        else:
            cls._python_constants["mu0"] = value

    # Gyromagnetic ratio
    @classmethod
    def gyromagnetic_ratio(cls) -> float:
        """Get gyromagnetic ratio [m/As]."""
        if cls._use_cpp():
            return _cpp_constants.PhysicalConstants.gyromagnetic_ratio()
        return cls._python_constants["gyromagnetic_ratio"]

    @classmethod
    def set_gyromagnetic_ratio(cls, value: float) -> None:
        """Set gyromagnetic ratio [m/As]."""
        if cls._use_cpp():
            _cpp_constants.PhysicalConstants.set_gyromagnetic_ratio(value)
        else:
            cls._python_constants["gyromagnetic_ratio"] = value

    # Gamma (frequency)
    @classmethod
    def gamma(cls) -> float:
        """Get gamma factor [Hz/T]."""
        return cls._python_constants["gamma"]

    @classmethod
    def set_gamma(cls, value: float) -> None:
        """Set gamma factor [Hz/T]."""
        if cls._use_cpp():
            _cpp_constants.PhysicalConstants.set_gyromagnetic_ratio(value)
        cls._python_constants["gamma"] = value

    # Gamma (angular frequency)
    @classmethod
    def gamma_rad(cls) -> float:
        """Get gamma factor [rad/(s*T)]."""
        return cls._python_constants["gamma_rad"]

    @classmethod
    def set_gamma_rad(cls, value: float) -> None:
        """Set gamma factor [rad/(s*T)].
        In CPP we use m/As, so we need to convert to rad/(s*T)
        by multiplying by 2*pi
        """
        cls._python_constants["gamma_rad"] = value

    # Tesla to A/m conversion
    @classmethod
    def TtoAm(cls) -> float:
        """Get Tesla to A/m conversion factor."""
        if cls._use_cpp():
            return _cpp_constants.PhysicalConstants.TtoAm()
        return cls._python_constants["TtoAm"]

    @classmethod
    def set_TtoAm(cls, value: float) -> None:
        """Set Tesla to A/m conversion factor."""
        if cls._use_cpp():
            _cpp_constants.PhysicalConstants.set_TtoAm(value)
        else:
            cls._python_constants["TtoAm"] = value

    @classmethod
    def AmtoT(cls) -> float:
        """Get A/m to Tesla conversion factor."""
        return 1.0 / cls.TtoAm()

    # Oersted conversions
    @classmethod
    def OetoAm(cls) -> float:
        """Get Oersted to A/m conversion factor."""
        if cls._use_cpp():
            return cls._python_constants["OetoAm"]  # Not in C++ constants yet
        return cls._python_constants["OetoAm"]

    @classmethod
    def set_OetoAm(cls, value: float) -> None:
        """Set Oersted to A/m conversion factor."""
        cls._python_constants["OetoAm"] = value
        cls._python_constants["AmtoOe"] = 1.0 / value

    @classmethod
    def AmtoOe(cls) -> float:
        """Get A/m to Oersted conversion factor."""
        return 1.0 / cls.OetoAm()

    # Planck's constant
    @classmethod
    def hplanck(cls) -> float:
        """Get Planck's constant [J⋅s]."""
        if cls._use_cpp():
            return _cpp_constants.PhysicalConstants.hbar() * 2 * math.pi
        return cls._python_constants["hplanck"]

    @classmethod
    def set_hplanck(cls, value: float) -> None:
        """Set Planck's constant [J⋅s]."""
        if cls._use_cpp():
            _cpp_constants.PhysicalConstants.set_hbar(value / (2 * math.pi))
        else:
            cls._python_constants["hplanck"] = value
            cls._python_constants["hbar"] = value / (2 * math.pi)

    @classmethod
    def hbar(cls) -> float:
        """Get reduced Planck's constant [J⋅s]."""
        if cls._use_cpp():
            return _cpp_constants.PhysicalConstants.hbar()
        return cls._python_constants["hbar"]

    @classmethod
    def set_hbar(cls, value: float) -> None:
        """Set reduced Planck's constant [J⋅s]."""
        if cls._use_cpp():
            _cpp_constants.PhysicalConstants.set_hbar(value)
        else:
            cls._python_constants["hbar"] = value
            cls._python_constants["hplanck"] = value * 2 * math.pi

    # Elementary charge
    @classmethod
    def echarge(cls) -> float:
        """Get elementary charge [C]."""
        if cls._use_cpp():
            return -_cpp_constants.PhysicalConstants.elementary_charge()  # Note: negative
        return cls._python_constants["echarge"]

    @classmethod
    def set_echarge(cls, value: float) -> None:
        """Set elementary charge [C]."""
        if cls._use_cpp():
            _cpp_constants.PhysicalConstants.set_elementary_charge(-value)  # Store positive in C++
        else:
            cls._python_constants["echarge"] = value

    # Electron mass
    @classmethod
    def me(cls) -> float:
        """Get electron mass [kg]."""
        return cls._python_constants["me"]  # Not in C++ constants yet

    @classmethod
    def set_me(cls, value: float) -> None:
        """Set electron mass [kg]."""
        cls._python_constants["me"] = value

    # Boltzmann constant
    @classmethod
    def boltzmann(cls) -> float:
        """Get Boltzmann constant [J/K]."""
        if cls._use_cpp():
            return _cpp_constants.PhysicalConstants.boltzmann_constant()
        return cls._python_constants["boltzmann"]

    @classmethod
    def set_boltzmann(cls, value: float) -> None:
        """Set Boltzmann constant [J/K]."""
        if cls._use_cpp():
            _cpp_constants.PhysicalConstants.set_boltzmann_constant(value)
        else:
            cls._python_constants["boltzmann"] = value

    # Derived constants
    @classmethod
    def bohr_magneton(cls) -> float:
        """Get Bohr magneton [J/T]."""
        return cls.echarge() * cls.hbar() / (2 * cls.me())

    # Utility methods
    @classmethod
    def reset_to_defaults(cls) -> None:
        """Reset all constants to their default values."""
        cls._python_constants = cls._FALLBACK_CONSTANTS.copy()
        if cls._use_cpp():
            _cpp_constants.PhysicalConstants.resetToDefaults()

    @classmethod
    def get_all_constants(cls) -> dict[str, float]:
        """Get a dictionary of all current constant values."""
        return {
            "mu0": cls.mu0(),
            "gyromagnetic_ratio": cls.gyromagnetic_ratio(),
            "gamma": cls.gamma(),
            "gamma_rad": cls.gamma_rad(),
            "TtoAm": cls.TtoAm(),
            "AmtoT": cls.AmtoT(),
            "OetoAm": cls.OetoAm(),
            "AmtoOe": cls.AmtoOe(),
            "hplanck": cls.hplanck(),
            "hbar": cls.hbar(),
            "echarge": cls.echarge(),
            "me": cls.me(),
            "boltzmann": cls.boltzmann(),
            "bohr_magneton": cls.bohr_magneton(),
        }


# For backward compatibility, provide module-level variables that dynamically
# get their values from the Constants class
def __getattr__(name: str):
    """Dynamic attribute access for backward compatibility."""
    if name == "mu0":
        return Constants.mu0()
    elif name == "gyromagnetic_ratio":
        return Constants.gyromagnetic_ratio()
    elif name == "gamma":
        return Constants.gamma()
    elif name == "gamma_rad":
        return Constants.gamma_rad()
    elif name == "TtoAm":
        return Constants.TtoAm()
    elif name == "AmtoT":
        return Constants.AmtoT()
    elif name == "OetoAm":
        return Constants.OetoAm()
    elif name == "AmtoOe":
        return Constants.AmtoOe()
    elif name == "hplanck":
        return Constants.hplanck()
    elif name == "hbar":
        return Constants.hbar()
    elif name == "echarge":
        return Constants.echarge()
    elif name == "me":
        return Constants.me()
    elif name == "boltzmann":
        return Constants.boltzmann()
    elif name == "bohr_magneton":
        return Constants.bohr_magneton()
    else:
        raise AttributeError(f"module '{__name__}' has no attribute '{name}'")


gyromagnetic_ratio = Constants.gyromagnetic_ratio()
gamma = Constants.gamma()
gamma_rad = Constants.gamma_rad()
TtoAm = Constants.TtoAm()
AmtoT = Constants.AmtoT()
OetoAm = Constants.OetoAm()
AmtoOe = Constants.AmtoOe()
hplanck = Constants.hplanck()
hbar = Constants.hbar()
echarge = Constants.echarge()
me = Constants.me()
boltzmann = Constants.boltzmann()
bohr_magneton = Constants.bohr_magneton()
mu0 = Constants.mu0()
# Export the Constants class and backward compatibility names
__all__ = [
    "Constants",
    "mu0",
    "gyromagnetic_ratio",
    "gamma",
    "gamma_rad",
    "TtoAm",
    "AmtoT",
    "OetoAm",
    "AmtoOe",
    "hplanck",
    "hbar",
    "echarge",
    "me",
    "boltzmann",
    "bohr_magneton",
]
