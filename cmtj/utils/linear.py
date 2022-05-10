from typing import Tuple

import numpy as np

from cmtj import CVector


class FieldScan:

    @staticmethod
    def _trig_compute(theta, phi) -> Tuple:
        st = np.sin(np.deg2rad(theta))
        ct = np.cos(np.deg2rad(theta))
        sp = np.sin(np.deg2rad(phi))
        cp = np.cos(np.deg2rad(phi))
        return st, ct, sp, cp

    @staticmethod
    def angle2vector(theta, phi, amplitude=1) -> CVector:
        st, ct, sp, cp = FieldScan._trig_compute(theta, phi)
        return CVector(
            st * cp * amplitude,
            st * sp * amplitude,
            ct * amplitude,
        )

    @staticmethod
    def vector2angle(x, y, z) -> Tuple:
        """
        https://github.com/numpy/numpy/issues/5228
        :returns (theta, phi, r)
        """
        r = np.sqrt(x**2 + y**2 + z**2)
        theta = np.rad2deg(np.arctan2(np.sqrt(x**2 + y**2), z))
        phi = np.rad2deg(np.arctan2(y, x))
        return theta, phi, r

    @staticmethod
    def cvector2angle(vector: CVector) -> Tuple:
        """
        https://github.com/numpy/numpy/issues/5228
        :returns (theta, phi, r)
        """
        return FieldScan.vector2angle(vector.x, vector.y, vector.z)

    @staticmethod
    def amplitude_scan(
        start: float,
        stop: float,
        steps: int,
        theta: float,
        phi: float,
        back: bool = False,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute a linear magnitude sweep. Angles given in deg.
        :param start:
        :param stop:
        :param steps:
        :param theta:
        :param phi:

        :returns: linear amplitude, field vectors
        """
        Hspan = np.linspace(start, stop, endpoint=True, num=steps)
        st, ct, sp, cp = FieldScan._trig_compute(theta, phi)
        Hx = st * cp * Hspan
        Hy = st * sp * Hspan
        Hz = ct * Hspan
        if back:
            forward = np.vstack((Hx, Hy, Hz)).T
            back = forward[::-1]
            return np.concatenate((Hspan, Hspan[::-1]),
                                  axis=0), np.concatenate((forward, back),
                                                          axis=0)
        return Hspan, np.vstack((Hx, Hy, Hz)).T

    @staticmethod
    def theta_scan(start: float, stop: float, steps: int, amplitude: float,
                   phi: float) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute a linear theta angle sweep. Angles given in deg.
        :param start:
        :param stop:
        :param steps:
        :param magnitude: magnitude of the scanned field.
        :param phi:
        """
        theta_span = np.linspace(start, stop, endpoint=True, num=steps)
        st, ct, sp, cp = FieldScan._trig_compute(theta_span, phi)
        Hx = st * cp * amplitude
        Hy = st * sp * amplitude
        Hz = ct * amplitude
        return theta_span, np.vstack((Hx, Hy, Hz)).T

    @staticmethod
    def phi_scan(start: float, stop: float, steps: int, amplitude: float,
                 theta: float) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute a linear phi angle sweep. Angles given in deg.
        :param start:
        :param stop:
        :param steps:
        :param magnitude: magnitude of the scanned field
        :param theta:
        """
        phi_span = np.linspace(start, stop, endpoint=True, num=steps)
        st, ct, sp, cp = FieldScan._trig_compute(theta, phi_span)
        Hx = st * cp * amplitude
        Hy = st * sp * amplitude
        Hz = ct * amplitude * np.ones_like(Hy)
        return phi_span, np.vstack((Hx, Hy, Hz)).T
