import numpy as np
from typing import Tuple


class FieldScan:
    @staticmethod
    def _trig_compute(theta, phi) -> Tuple:
        st = np.sin(np.deg2rad(theta))
        ct = np.cos(np.deg2rad(theta))
        sp = np.sin(np.deg2rad(phi))
        cp = np.cos(np.deg2rad(phi))
        return st, ct, sp, cp

    @staticmethod
    def amplitude_scan(start: float, stop: float, steps: int, theta: float,
                       phi: float) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute a linear magnitude sweep.
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
        return Hspan, np.vstack((Hx, Hy, Hz)).T

    @staticmethod
    def theta_scan(start: float, stop: float, steps: int, amplitude: float,
                   phi: float) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute a linear theta angle sweep.  
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
        Compute a linear phi angle sweep.  
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
