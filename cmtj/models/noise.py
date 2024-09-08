from typing import Literal

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm


def noise_model(
    N: int,
    steps: int = 3e5,
    thermal_noise_std: float = 1e-3,
    background_thermal_noise_std: float = 1e-3,
    amplitude: float = 1e-3,
    enable_oscillations: bool = False,
    volume_distribution: Literal["pareto", "uniform"] = "pareto",
    volume_distribution_params: dict = None,
    freq_distribution: Literal["uniform", "functional"] = "uniform",
    freq_distribution_params: dict = None,
    frequency_scale: float = 1000,
    time_scale: float = 1e-9,
    phase_std: float = np.pi / 12,
    dims: int = 1,
    seed: int = 42,
    offset: int = 1,
    save_vectors: bool = False,
    verbose: bool = False,
    N_background_scale: int = 10,
):
    """
    Generate a basic noise model.

    :param N: Number of elements in the noise model.
    :type N: int
    :param steps: Number of simulation steps, defaults to 3e5.
    :type steps: int, optional
    :param thermal_noise_std: Standard deviation of the thermal noise, defaults to 1e-3.
    :type thermal_noise_std: float, optional
    :param amplitude: Amplitude of the oscillating frequencies, defaults to 1e-3.
    :type amplitude: float, optional
    :param enable_oscillations: Enable oscillations, defaults to False.
    :type enable_oscillations: bool, optional
    :param volume_distribution: Volume distribution type, defaults to "pareto".
    :type volume_distribution: str, optional
    :param volume_distribution_params: Parameters for the volume distribution, defaults to None.
    :type volume_distribution_params: dict, optional
    :param freq_distribution: Frequency distribution type, defaults to "uniform".
    :type freq_distribution: str, optional
    :param freq_distribution_params: Parameters for the frequency distribution, defaults to None.
    :type freq_distribution_params: dict, optional
    :param frequency_scale: Frequency scale, defaults to 1000.
    :type frequency_scale: float, optional
    :param time_scale: Time scale, defaults to 1e-9.
    :type time_scale: float, optional
    :param phase_std: Standard deviation of the phase, defaults to np.pi / 12.
    :type phase_std: float, optional
    :param dims: Number of dimensions, defaults to 1.
    :type dims: int, optional
    :param seed: Random seed, defaults to 42.
    :type seed: int, optional
    :param offset: Offset for the simulation, defaults to 1.
    :type offset: int, optional
    :param save_vectors: Save the stepwise m vectors per each volume, defaults to False.
    :type save_vectors: bool, optional
    :param verbose: Enable verbose mode, defaults to False.
    :type verbose: bool, optional
    :param N_background_scale: Number of background domains for oscillations and
        background noise, defaults to 10xN.
    :type N_background_scale: int, optional
    :return: A tuple containing the following elements:
        - m_values (ndarray): Array of shape (steps, dims) representing the simulated values.
        - volumes (ndarray): Array of shape (N, 1) representing the volumes.
        - freqs (ndarray): Array of shape (N,) representing the frequencies.
        - time_scale (float): The time scale used in the simulation.
    :rtype: tuple
    """
    if volume_distribution_params is None:
        volume_distribution_params = {"shape": 1, "scale": 0.05}
    if freq_distribution_params is None:
        freq_distribution_params = {"low": 1, "high": 5000}
    rng = np.random.default_rng(seed=seed)
    if volume_distribution == "pareto":
        volumes = rng.gamma(**volume_distribution_params, size=N)
    elif volume_distribution == "uniform":
        volumes = rng.random(N)
        volumes = volumes / np.sum(volumes)

    steps = int(steps)
    volumes = volumes[:, np.newaxis]
    volumes = np.sort(volumes, axis=0)

    if freq_distribution == "functional":
        freqs = (frequency_scale / volumes.ravel()).astype(int)[::-1] + 1
    elif freq_distribution == "uniform":
        freqs = rng.integers(**freq_distribution_params, size=N)

    Np = N_background_scale * N
    freqs_osc = rng.uniform(10, 1e11, Np)
    phases = np.zeros(Np)
    if phase_std > 0:
        phases = rng.random(Np) * phase_std

    vector_values = rng.random((N, dims))
    m_values = np.zeros((steps, dims))
    f_counts = np.zeros_like(freqs)
    vectors = []
    triggers = 0

    def _oscillations(i: int):
        return amplitude * np.sin(2 * np.pi * freqs_osc * i * time_scale + phases).reshape(-1, 1)

    def _background_noise(i: int):
        return rng.normal(0, background_thermal_noise_std, dims)

    if enable_oscillations and background_thermal_noise_std > 0:
        raise ValueError(
            "Cannot have both oscillations and background thermal noise enabled."
            " Either set enable oscillations = False or background thermal noise = 0."
        )
    if enable_oscillations:
        osc_fn = _oscillations
    elif background_thermal_noise_std > 0:
        osc_fn = _background_noise
    for i in tqdm(range(offset, steps + offset), total=steps):
        osc_vals = osc_fn(i).sum()
        freq_disturbance = rng.integers(-3, 3, N)
        freq_mask = i % (freqs + freq_disturbance) == 0
        fsum = np.sum(freq_mask)
        f_counts[freq_mask] += 1
        if fsum > 0:
            triggers += 1
            vector_values[freq_mask] = rng.normal(0, thermal_noise_std, (fsum, dims))
            m_values[i - offset] += np.sum(volumes * vector_values, axis=0) + osc_vals
        else:
            m_values[i - offset] = osc_vals + m_values[i - offset - 1]

        if save_vectors:
            vectors.append(vector_values.copy())

    if verbose:
        print(f"Triggers {triggers} out of {steps} steps")
    if save_vectors:
        return m_values, volumes, freqs, time_scale, f_counts, vectors
    return m_values, volumes, freqs, time_scale, f_counts


def autocorrelation(x, dT):
    """
    Compute the autocorrelation of the signal, based on the properties of the
    power spectral density of the signal.
    Taken from the StackOverflow answer:
    https://stackoverflow.com/questions/643699/how-can-i-use-numpy-correlate-to-do-autocorrelation
    :param x: the signal
    :param dT: the time step
    """
    xp = x - np.mean(x)
    f = np.fft.fft(xp)
    p = np.abs(f) ** 2
    pi = np.fft.ifft(p)
    autocorr = np.real(pi)[: x.size // 2] / np.sum(xp**2)

    # Create a lag array
    lag = np.arange(0, len(autocorr)) * dT

    return lag, autocorr


def plot_noise_data(m_values: np.ndarray, volumes: np.ndarray, freqs: np.ndarray, time_scale: float):
    """
    Plot noise data:
    - Autocorrelation
    - Volume vs Frequency
    - Frequency vs Power
    - Time vs Amplitude

    :param m_values: An array of magnetic values.
    :type m_values: np.ndarray
    :param volumes: An array of volumes.
    :type volumes: np.ndarray
    :param freqs: An array of frequencies.
    :type freqs: np.ndarray
    :param time_scale: The time scale.
    :type time_scale: float

    :returns: The generated figure.
    :rtype: matplotlib.figure.Figure
    """
    dims = m_values.shape[1]
    R = m_values.dot(np.array([1, 0, 0])) if dims == 3 else m_values
    lag, autocorr = autocorrelation(R.ravel(), time_scale)
    k = m_values.shape[0]
    with plt.style.context(["science", "nature"]):
        w, h = plt.figaspect(1.0 / 4.0)
        fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, dpi=300, figsize=(w, h))
        ax1.plot(lag, autocorr, color="royalblue")
        ax1.set_xlabel("Time lag (s)")
        ax1.set_ylabel("Autocorrelation")
        ax2.plot(volumes, freqs / 1000, color="crimson")
        # histogram of volumes
        ax25 = ax2.twinx()
        ax25.hist(volumes, bins=min(100, len(volumes)), color="navy", alpha=0.5, label="Count")
        ax25.set_ylabel("Count", rotation=-90, labelpad=10)
        ax25.legend()
        ax2.set_xlabel("Area (a.u.)")
        ax2.set_ylabel("Modulo step activation (1000x)")
        y = np.fft.fft(m_values, axis=0)
        y = np.power(np.abs(y), 2)
        y = y[: int(k // 2)]
        x = np.fft.fftfreq(int(k), time_scale)
        x = x[: int(k // 2)]
        ax3.plot(x, y, color="royalblue")
        ax3.set_xscale("log")
        ax3.set_yscale("log")
        ax3.set_xlabel("Frequency (Hz)")
        ax3.set_ylabel("Power (a.u.)")
        x_base = np.arange(0, k * time_scale, time_scale)
        ax4.plot(x_base, m_values, color="forestgreen")
        ax4.set_xlabel("Time")
        ax4.set_ylabel("Amplitude")
        fig.subplots_adjust(wspace=0.55)

        # add letters
        import matplotlib.transforms as mtransforms

        for label, ax in zip("abcd", (ax1, ax2, ax3, ax4)):
            # label physical distance in and down:
            trans = mtransforms.ScaledTranslation(10 / 72, -5 / 72, fig.dpi_scale_trans)
            ax.text(
                0.0,
                1.0,
                f"{label})",
                transform=ax.transAxes + trans,
                # fontsize="medium",
                verticalalignment="top",
                color="black",
                bbox=dict(facecolor="none", alpha=0.4, edgecolor="none", pad=3.0),
            )

    return fig


def create_noise_animation(
    volumes: list,
    vector_values: np.ndarray,
    save_type: str = "gif",
    max_frames: int = 1000,
):
    """
    Create a 2D noise animation given volumes/areas and vector values.

    :param volumes: A list of volumes.
    :type volumes: list
    :param vector_values: An array of vector values.
    :type vector_values: np.ndarray
    :param save_type: The type of animation to save. Defaults to "gif".
    :type save_type: str, optional
    :param max_frames: The maximum number of frames. Defaults to 1000.
    :type max_frames: int, optional
    :returns: None
    """

    rng = np.random.default_rng(seed=42)
    vector_values = np.asarray(vector_values).squeeze()
    v = volumes.ravel()
    v = v / v.sum()
    n = 1000
    xx, yy = np.meshgrid(np.arange(0, n), np.arange(0, n))
    values = np.zeros_like(xx, dtype=np.float32)
    volume_masks = []
    for i, volume in enumerate(v):
        x0, y0 = rng.integers(0, n, 2)
        shape = (xx - x0) ** 2 + (yy - y0) ** 2
        mask = shape <= (volume / np.pi) * ((n / 2) ** 2)
        values[mask] = vector_values[105, i]
        volume_masks.append(mask)

    # create matplotlib animation
    fig, ax = plt.subplots(dpi=300)
    im = ax.imshow(values, cmap="viridis", vmin=-1, vmax=1)
    # remove axes
    ax.axis("off")

    def update(step):
        for i, mask in enumerate(volume_masks):
            values[mask] = vector_values[step, i]
        im.set_array(values)
        ax.set_title(f"Step {step}")

    ani = animation.FuncAnimation(fig, update, frames=max_frames, interval=0.8)
    if save_type == "gif":
        ani.save("animation.gif", writer="pillow")
    else:
        ani.save("animation.mp4", writer="ffmpeg")
