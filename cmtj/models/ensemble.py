import numpy as np


def meinert_model(phi, V1, V2, phase_offset, offset):
    """
    Fits to Meinert model.
    :param phi: angle in degrees, given parameter
    :param V1: (Hfl/Hext), fitting parameter
    :param V2: (Va*Hdl/Heff), fitting parameter
    :param phase_offset: phase offset in degrees, fitting parameter
    :param offset: offset in V, fitting parameter
    V2omega = (Acos(2(phi - phase_offset)) - B)*cos(phi - phase_offset) + offset
    """
    deg_rad = np.deg2rad(phi - phase_offset)
    return (V1 * np.cos(2 * (deg_rad)) - V2) * np.cos(deg_rad) + offset


def symmetric_lorentz(H, dH, Hr, Vs):
    """
    Symmetric Lorentzian function.
    :param H: applied field in A/m
    :param dH: half width at half maximum in A/m
    :param Hr: resonance field in A/m
    """
    dH2 = dH**2
    return Vs * dH2 / ((H - Hr) ** 2 + dH2)


def antisymmetric_lorentz(H, dH, Hr, Vas):
    """
    Antisymmetric Lorentzian function.
    :param H: applied field in A/m
    :param dH: half width at half maximum in A/m
    :param Hr: resonance field in A/m
    """
    dH2 = dH**2
    dHr = H - Hr
    return Vas * dH * dHr / (np.power(dHr, 2) + dH2)


def mixed_lorentz(H, dH, Hr, Va, Vas):
    """
    Mixed Lorentzian function.
    :param H: applied field in A/m
    :param dH: half width at half maximum in A/m
    :param Hr: resonance field in A/m
    :param Va: amplitude of symmetric Lorentzian
    :param Vas: amplitude of antisymmetric Lorentzian
    """
    return symmetric_lorentz(H, dH, Hr, Va) + antisymmetric_lorentz(H, dH, Hr, Vas)
