import numpy as np
from scipy.optimize import root


class RootFinder:
    """Adopted from: https://stackoverflow.com/a/65185377/3588442"""

    def __init__(self, start, stop, step=0.01, root_dtype="float32", xtol=1e-9):
        self.start = start
        self.stop = stop
        self.step = step
        self.xtol = xtol
        self.roots = np.array([], dtype=root_dtype)

    def add_to_roots(self, x):
        if (x < self.start) or (x > self.stop):
            return  # outside range
        if any(abs(self.roots - x) < self.xtol):
            return  # root already found.

        self.roots = np.append(self.roots, x)

    def find(self, f, *args, **kwargs):
        current = self.start
        nsteps = int((self.stop - self.start) / self.step)
        for x0 in np.linspace(self.start, self.stop + self.step, nsteps):
            if x0 < current:
                continue
            x = self.find_root(f, x0, *args, **kwargs)
            if x is None:  # no root found.
                continue
            current = x
            self.add_to_roots(x)

        return self.roots

    def find_root(self, f, x0, *args, fprime=None):
        sol = root(f, x0, args=args, jac=fprime, method="hybr")
        return sol.x[0] if sol.success else None
