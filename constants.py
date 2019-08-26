import numpy as np


class Constants:
    @property
    def MAGNETIC_PERMEABILITY(self):
        return 12.57e-7

    # @property
    # def GYRO(self):
    #     return 176.08596e9
    @property
    def GYRO(self):
        return 2.21e5

    @property
    def DAMPING(self):
        return 0.02

    @property
    def TtoAm(self):
        return 795774.715459

    @property
    def HBAR(self):
        return 6.62607015e-34/(2*np.pi)

    @property
    def ELECTRON_CHARGE(self):
        return 1.60217662e-19
