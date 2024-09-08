import math
from dataclasses import dataclass
from typing import Union

import numpy as np

from cmtj import CVector


@dataclass
class VectorObj:
    """Vector object for standard manipulation.
    Alternative to CVectors (which are used in the C++ code).
    Easier to modify and manipulate, but slower.
    :param theta: positive z-axis angle (in xz plane) in radians.
    :param phi: positive x-axis (in xy plane) angle in radians
    :param mag: magnitude of the vector, if not set defaults to 1 *unit vector*
    """

    theta: float  # in radians
    phi: float  # rad
    mag: float = 1

    def __add__(self, other):
        """Adds two vectors"""
        return VectorObj.from_cvector(self.to_cvector() + other.to_cvector())

    def __mul__(self, other: Union["VectorObj", float]):
        """Multiplies a vector by a scalar"""
        if isinstance(other, VectorObj):
            return self._componentwise_mul(other)
        return VectorObj.from_cvector(self.to_cvector() * other)

    def __rmul__(self, other: Union["VectorObj", float]):
        """Multiplies a vector by a scalar"""
        return self.__mul__(other)

    def __repr__(self) -> str:
        return f"VectorObj(theta={self.theta}, phi={self.phi}, mag={self.mag})"

    def __hash__(self) -> int:
        return hash(str(self))

    def __eq__(self, __value: "VectorObj") -> bool:
        return self.theta == __value.theta and self.phi == __value.phi and self.mag == __value.mag

    def _componentwise_mul(self, other):
        coors = self.get_cartesian()
        other_coords = other.get_cartesian()
        return VectorObj.from_cartesian(
            coors[0] * other_coords[0],
            coors[1] * other_coords[1],
            coors[2] * other_coords[2],
        )

    def get_cartesian(self):
        """Returns the vector in Cartesian coordinates with (x, y, z) compnents"""
        return VectorObj.from_spherical(self.theta, self.phi, self.mag)

    @staticmethod
    def from_spherical(theta, phi, mag=1):
        """Creates a Cartesian vector from spherical components"""
        return [
            mag * math.sin(theta) * math.cos(phi),
            mag * math.sin(theta) * math.sin(phi),
            mag * math.cos(theta),
        ]

    @staticmethod
    def from_cartesian(x: float, y: float, z: float):
        """Creates a spherical vector from Cartesian components"""
        mag = math.sqrt(x**2 + y**2 + z**2)
        if mag == 0:
            return VectorObj(0, 0, 0)
        theta = math.acos(z / mag)
        phi = math.atan2(y, x)
        return VectorObj(theta, phi, mag)

    @staticmethod
    def from_cvector(cvector: CVector):
        """Creates a spherical vector from Cartesian components"""
        mag = cvector.length()
        if mag == 0:
            return VectorObj(0, 0, 0)
        theta = math.acos(cvector.z / mag)
        phi = math.atan2(cvector.y, cvector.x)
        return VectorObj(theta, phi, mag)

    def to_cvector(self):
        """Creates a Cartesian vector from spherical components"""
        return CVector(*self.get_cartesian())


def box_muller_random(mean, std):
    """
    Generates Gaussian noise with mean and standard deviation
    using the Box-Muller transform.
    https://en.wikipedia.org/wiki/Box–Muller_transform
    :param mean: mean of the Gaussian.
    :param std: standard deviation of the Gaussian.
    """
    u1 = np.random.uniform(0, 1)
    u2 = np.random.uniform(0, 1)
    mag = std * math.sqrt(-2.0 * math.log(u1))
    z0 = mag * math.cos(2 * math.pi * u2) + mean
    z1 = mag * math.sin(2 * math.pi * u2) + mean
    return z0, z1


def perturb_position(eq_point, pmax=1e-3):
    """
    Perturbs an equilibrium point by a random amount.
    """
    return np.asarray(eq_point) + np.random.normal(0, pmax, len(eq_point))
