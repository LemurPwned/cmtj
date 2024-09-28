import math

from numba import njit


def constant_driver(t, amp):
    return amp


@njit
def sine_driver(t, amp, freq, phase):
    return amp * math.sin(freq * t + phase)


@njit
def decay_driver(t, amp, tau):
    return amp * math.exp(-t / tau)


@njit
def gaussian_driver(t, amp, sigma):
    return amp * math.exp(-(t**2) / sigma**2)
