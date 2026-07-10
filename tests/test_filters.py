import numpy as np

from cmtj.utils.filters import Filters


def test_bandpass_filter_accepts_hertz_frequencies():
    data = np.sin(np.linspace(0, 2 * np.pi, 1000))

    filtered = Filters.butter_bandpass_filter(data, (100, 200), fs=1000)

    assert filtered.shape == data.shape
    assert np.isfinite(filtered).all()
