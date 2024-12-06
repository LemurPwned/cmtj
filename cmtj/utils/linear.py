import numpy as np

from cmtj import CVector


class FieldScan:
    @staticmethod
    def _trig_compute(theta, phi) -> tuple:
        """Compute trigonometric functions for theta and phi.
        :param theta: theta angle in [deg].
        :param phi: phi angle in [deg].
        :returns: trigonometric functions for theta and phi."""
        st = np.sin(np.deg2rad(theta))
        ct = np.cos(np.deg2rad(theta))
        sp = np.sin(np.deg2rad(phi))
        cp = np.cos(np.deg2rad(phi))
        return st, ct, sp, cp

    @staticmethod
    def angle2vector(theta, phi, amplitude=1) -> CVector:
        """Convert spherical coordinates to cartesian coordinates.
        :param theta: polar angle in degrees.
        :param phi: azimuthal angle in degrees.
        :param amplitude: amplitude of target vector.
        :returns: cartesian vector."""
        st, ct, sp, cp = FieldScan._trig_compute(theta, phi)
        return CVector(
            st * cp * amplitude,
            st * sp * amplitude,
            ct * amplitude,
        )

    @staticmethod
    def vector2angle(x, y, z) -> tuple:
        """Convert cartesian coordinates to spherical coordinates.
        :param x: x coordinate of the vector.
        :param y: y coordinate of the vector.
        :param z: z coordinate of the vector.
        :returns (theta, phi, r)
        https://github.com/numpy/numpy/issues/5228
        """
        r = np.sqrt(x**2 + y**2 + z**2)
        theta = np.rad2deg(np.arctan2(np.sqrt(x**2 + y**2), z))
        phi = np.rad2deg(np.arctan2(y, x))
        return theta, phi, r

    @staticmethod
    def cvector2angle(vector: CVector) -> tuple:
        """
        :param vector: cartesian vector.
        :returns (theta, phi, r)
        https://github.com/numpy/numpy/issues/5228
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
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Compute a linear magnitude sweep. Angles given in deg.
        :param start: start of the sweep
        :param stop: end of the sweep
        :param steps: number of steps
        :param theta: polar angle in deg.
        :param phi: azimuthal angle in deg.
        :returns: linear amplitude, field vectors
        """
        Hspan = np.linspace(start, stop, endpoint=True, num=steps)
        st, ct, sp, cp = FieldScan._trig_compute(theta, phi)
        Hx = st * cp * Hspan
        Hy = st * sp * Hspan
        Hz = ct * Hspan
        if back:
            forward = np.vstack((Hx, Hy, Hz)).T
            return np.concatenate((Hspan[:-1], Hspan[::-1]),
                                  axis=0), np.concatenate(
                                      (forward[:-1], forward[::-1]), axis=0)
        return Hspan, np.vstack((Hx, Hy, Hz)).T

    @staticmethod
    def theta_scan(
        start: float, stop: float, steps: int, amplitude: float, phi: float
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Compute a linear theta angle sweep. Angles given in deg.
        :param start: polar angle start of the sweep
        :param stop: polar angle end of the sweep
        :param steps: number of steps
        :param amplitude: amplitude of the scanned field.
        :param phi: azimuthal angle in deg.
        """
        theta_span = np.linspace(start, stop, endpoint=True, num=steps)
        st, ct, sp, cp = FieldScan._trig_compute(theta_span, phi)
        Hx = st * cp * amplitude
        Hy = st * sp * amplitude
        Hz = ct * amplitude
        return theta_span, np.vstack((Hx, Hy, Hz)).T

    @staticmethod
    def phi_scan(
        start: float, stop: float, steps: int, amplitude: float, theta: float
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Compute a linear phi angle sweep. Angles given in deg.
        :param start: azimuthal angle start of the sweep
        :param stop: azimuthal angle end of the sweep
        :param steps: number of steps
        :param amplitude: amplitude of the scanned field
        :param theta: polar angle in deg.
        """
        phi_span = np.linspace(start, stop, endpoint=True, num=steps)
        st, ct, sp, cp = FieldScan._trig_compute(theta, phi_span)
        Hx = st * cp * amplitude
        Hy = st * sp * amplitude
        Hz = ct * amplitude * np.ones_like(Hy)
        return phi_span, np.vstack((Hx, Hy, Hz)).T
