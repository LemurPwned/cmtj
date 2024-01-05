from typing import List

import numpy as np

GENERIC_BOUNDS = {
    "Ms": (0.2, 2.5),  # T
    "K": (0.1, 3e3),  # kJ/m^3
    "J": (-5e3, 5e3),  # uJ/m^2
}

GENERIC_UNITS = {
    "Ms": "T",
    "K": "kJ/m^3",
    "J": "uJ/m^2",
}

def extract_max_resonance_lines(
    spectrum,
    h_vals: List[float],
    frequencies: List[float],
    N_layers: int,
    mask_size_x=0,
    mask_size_y=0,
):
    masked_spectrum = spectrum[mask_size_x:, mask_size_y:]
    frequencies = frequencies[mask_size_y:].squeeze()
    h_vals = h_vals[mask_size_x:]
    linespectra = {}
    for i in range(len(h_vals)):
        grad = np.gradient(masked_spectrum[i, :])
        sign = np.sign(grad)
        diff = np.diff(sign)
        indices = np.where(diff != 0)[0]
        grad = grad / np.max(np.abs(grad))

        # take two largest values at the indices
        amp_vals = masked_spectrum[i, indices]
        amp_max_indx = np.argsort(amp_vals)
        indices = indices[amp_max_indx]

        indices = indices[-N_layers:]
        n_diff = N_layers - len(indices)
        freqs = frequencies[indices].ravel() + [None] * n_diff
        linespectra[h_vals[i]] = freqs
    return linespectra
