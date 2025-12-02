"""
Probabilistic Bit (p-bit) Stochastic Magnetization Dynamics Simulation

GOAL:
This simulation studies probabilistic computing elements (p-bits) based on magnetic
tunnel junctions operating in the superparamagnetic regime. The simulation analyzes
stochastic magnetization switching under thermal fluctuations and spin-transfer torque,
characterizing the probabilistic behavior essential for neuromorphic computing applications.

CRITICAL SIMULATION MECHANISMS:
1. Stochastic Dynamics: Implements thermal fluctuations using Heun solver with finite
   temperature to capture random magnetization switching events
2. Spin-Transfer Torque (STT): Models current-induced torque with Slonczewski formalism
   including both parallel and perpendicular torque components
3. Uses reduced magnetic volume and moderate anisotropy to
   enable thermal activation over energy barriers
4. Statistical Sampling: Performs multiple simulation runs with different initial
   conditions to build statistical distributions of switching behavior
5. Binary State Analysis: Converts continuous magnetization trajectories to binary
   states using hysteresis thresholding for digital logic applications
6. Autocorrelation Analysis: Computes temporal correlations to characterize switching
   rates and memory effects in the stochastic dynamics
7. Field-Dependent Statistics: Maps switching probability and dwell times as functions
   of applied magnetic field and current bias

The simulation generates probabilistic switching statistics, autocorrelation functions,
and field-dependent probability maps that characterize the p-bit performance. This
enables optimization of device parameters for neuromorphic computing, where controlled
randomness and tunable switching probabilities are essential for implementing
probabilistic algorithms and neural network functions.
"""

import contextlib

import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import mpl_toolkits.mplot3d.art3d as art3d
import numpy as np
from matplotlib.patches import Circle
from tqdm import tqdm

from cmtj import (
    AxialDriver,
    CVector,
    Junction,
    Layer,
    SolverMode,
    constantDriver,
)

with contextlib.suppress(ImportError):
    import scienceplots  # noqa: F401

r = 1
pi = np.pi
cos = np.cos
sin = np.sin
phi, theta = np.mgrid[0.0:pi:100j, 0.0 : 2.0 * pi : 100j]
x = r * sin(phi) * cos(theta)
y = r * sin(phi) * sin(theta)
z = r * cos(phi)


def plot_trajectories(log: dict[str, list[float]], title: str):
    with plt.style.context(["science", "no-latex"]):
        fig = plt.figure(figsize=(12, 6))
        ax = fig.add_subplot(1, 2, 1, projection="3d")
        m = np.asarray([log["free_mx"], log["free_my"], log["free_mz"]])
        ax.plot3D(m[0], m[1], m[2], color="blue")
        ax.set_axis_off()
        ax.plot_surface(x, y, z, rstride=2, cstride=2, color="c", alpha=0.3, linewidth=0.1)
        ax.scatter([0], [0], [1], color="crimson", alpha=1.0, s=50)
        ax2 = fig.add_subplot(1, 2, 2)
        R = np.asarray(log["R"])
        ax2.plot(log["time"], R, alpha=0.3)
        hyst = R
        hyst[hyst < 150] = 100
        hyst[hyst >= 150] = 200
        ax2.plot(np.asarray(log["time"]) * 1e9, hyst, color="r")
        ax2.set_xlabel("Time [ns]")
        ax2.set_ylabel("Resistance [Ohm]")
        fig.suptitle(title)
        fig.tight_layout()
        return fig


def measure_event_time(series):
    pos_event_times, neg_event_times = [], []
    cumevent = 0
    for n in range(len(series) - 1):
        if series[n] == series[n + 1] == 1 or series[n] == series[n + 1] == 0:
            cumevent += 1
        elif series[n] == 1 and series[n + 1] == 0:
            # transition 1->0
            if cumevent > 1:  # do not count single dt switch
                pos_event_times.append(cumevent)
            cumevent = 0
        elif series[n] == 0 and series[n + 1] == 1:
            # transition 0->1
            if cumevent > 1:  # do not count single dt switch
                neg_event_times.append(cumevent)
            cumevent = 0
    return np.asarray(pos_event_times), np.asarray(neg_event_times)


def hyst(x, th_lo, th_hi, initial=False):
    """
    x : Numpy Array
        Series to apply hysteresis to.
    th_lo : float or int
        Below this threshold the value of hyst will be False (0).
    th_hi : float or int
        Above this threshold the value of hyst will be True (1).
    """

    if th_lo > th_hi:  # If thresholds are reversed, x must be reversed as well
        x = x[::-1]
        th_lo, th_hi = th_hi, th_lo
        rev = True
    else:
        rev = False

    hi = x >= th_hi
    lo_or_hi = (x <= th_lo) | hi

    ind = np.nonzero(lo_or_hi)[0]  # Index für alle darunter oder darüber
    if not ind.size:  # prevent index error if ind is empty
        x_hyst = np.zeros_like(x, dtype=bool) | initial
    else:
        cnt = np.cumsum(lo_or_hi)  # from 0 to len(x)
        x_hyst = np.where(cnt, hi[ind[cnt - 1]], initial)

    if rev:
        x_hyst = x_hyst[::-1]

    return x_hyst


def autocorrelation(x):
    """
    Compute the autocorrelation of the signal, based on the properties of the
    power spectral density of the signal.
    """
    xp = x - np.mean(x)
    f = np.fft.fft(xp)
    p = np.abs(f) ** 2
    pi = np.fft.ifft(p)
    return np.real(pi)[: x.size // 2] / np.sum(xp**2)


def discreticise_space(signal: np.ndarray):
    signal = np.asarray(signal)
    # normalised = (signal - signal.min()) / (signal.max() - signal.min())
    normalised = (signal - 100) / 100
    hysteresis = hyst(normalised, 0.4, 0.6)

    auto_correlation = autocorrelation(normalised)

    return normalised, hysteresis, auto_correlation


def transition_counter(signal: np.ndarray):
    h = signal.astype(int)
    rh = np.roll(signal, shift=1).astype(int)
    return np.abs(h - rh).sum()


def mean_one_time(time, binary_signal):
    one_state = np.argwhere(binary_signal == 1).ravel()
    return time[one_state].sum() / time.sum()


def sample_initial_conditions(theta0, phi0, eps):
    """Sample initial conditions for the stochastic simulation
    from starting point + radius"""
    theta0 = np.random.uniform(0, theta0 + eps)
    phi0 = np.random.uniform(0, 2 * np.pi)
    return CVector.fromSpherical(theta0, phi0, 1.0)


eps = 0.3
theta0 = 0.0
phi0 = 0.0
initial_conditions = [sample_initial_conditions(theta0, phi0, eps) for _ in range(100)]

demagTensor = [CVector(0.0, 0.0, 0.0), CVector(0.0, 0.0, 0.0), CVector(0.0, 0.0, 1.0)]

damping = 0.1
currentDensity = 1e10
beta = 0.9
spinPolarisation = 1.0

l1 = Layer.createSTTLayer(
    id="free",
    mag=CVector(1.0, 0.0, 0.0),
    anis=CVector(1, 0.0, 0.0),
    Ms=1.6,
    thickness=2.0e-9,
    cellSurface=7e-10 * 7e-10 * np.pi,
    demagTensor=demagTensor,
    damping=damping,
    SlonczewskiSpacerLayerParameter=1.0,
    spinPolarisation=spinPolarisation,
    beta=beta,
)

l1.setReferenceLayer(CVector(1, 0.0, 0.0))
junction = Junction([l1], 100, 200)

junction.setLayerAnisotropyDriver("free", constantDriver(350e3))
# current driver
junction.setLayerTemperatureDriver("free", constantDriver(10))
junction.setLayerCurrentDriver("free", constantDriver(-currentDensity))

field_trajs = []
field_vec = (-500e3, -50e3, 0, 50e3, 500e3)
for field in field_vec:
    trajectories = []
    junction.setLayerExternalFieldDriver(
        "free",
        AxialDriver(field, 0, 0),
    )

    for ic in tqdm(initial_conditions):
        junction.clearLog()
        junction.setLayerMagnetisation("free", ic)
        junction.runSimulation(5e-9, 1e-13, 1e-13, False, False, solverMode=SolverMode.Heun)
        log = junction.getLog()
        m = np.asarray([log["free_mx"], log["free_my"], log["free_mz"]])
        trajectories.append(m)
    field_trajs.append(np.asarray(trajectories))


with plt.style.context(["science", "nature"]):
    fig, ax = plt.subplots(
        1,
        len(field_vec),
        sharey=True,
        dpi=300,
        figsize=plt.figaspect(1.0 / len(field_vec)),
    )

    counts, used_bins = np.histogram(np.asarray(field_trajs[2])[:, 0, -1], bins=50)
    kspace = 1

    # color scale for 3D plot
    # scalar mappling
    norm = plt.Normalize(-1, 1)
    cmap = plt.cm.viridis

    for i in range(len(field_vec)):
        final_states = np.asarray(field_trajs[i])[:, 0, -1]
        # ax[i].hist(final_states, bins=used_bins, density=False, color="royalblue")

        hist, bins = np.histogram(final_states, bins=used_bins, density=False)
        bin_centers = 0.5 * (bins[:-1] + bins[1:])
        colors = cmap(norm(bin_centers))

        for bin_start, bin_end, color in zip(bins[:-1], bins[1:], colors):
            ax[i].bar(
                bin_start,
                hist[np.where(bins[:-1] == bin_start)[0][0]],
                width=bin_end - bin_start,
                color=color,
                edgecolor="black",
            )

        ax[i].set_xlim([-1.1, 1.1])

        # 3D plot inset
        inset = ax[i].inset_axes([-1.1, 0, 2.2, 100], transform=ax[i].transData, projection="3d")
        inset.set_axis_off()
        inset.patch.set_alpha(0.0)
        inset.plot_surface(x, y, z, rstride=1, cstride=1, color="azure", alpha=0.1, linewidth=0.1)
        ax[i].set_title(r"$\mathrm{H}_\mathrm{x} = $" + f"{field_vec[i] / 1e3} kA/m")
        for m in field_trajs[i][::kspace]:
            colors = cmap(norm(m[0][-1]))
            inset.scatter(
                m[0][-1],
                m[1][-1],
                m[2][-1],
                color=colors,
                edgecolor="none",
                alpha=1.0,
                s=6,
            )

            inset.scatter(
                m[0][0],
                m[1][0],
                m[2][0],
                color="crimson",
                edgecolor="black",
                alpha=0.8,
                s=5,
            )
        p = Circle(
            (0, 0),
            np.sin(eps),
            alpha=0.3,
            color="crimson",
        )
        inset.add_patch(p)
        art3d.pathpatch_2d_to_3d(p, z=np.cos(eps), zdir="z")

        ax[i].set_xlabel(r"$\mathrm{m}_x$")
    ax[0].set_ylabel(r"\# of trajectories")
    fig.subplots_adjust(wspace=0.05)
    for label, ax_ in zip("abcdefg", ax.flatten()):
        # label physical distance in and down:
        trans = mtransforms.ScaledTranslation(10 / 72, -5 / 72, fig.dpi_scale_trans)
        ax_.text(
            0.05,
            1.0,
            f"{label})",
            transform=ax_.transAxes + trans,
            fontsize="medium",
            verticalalignment="top",
            color="black",
            bbox=dict(facecolor="none", alpha=0.4, edgecolor="none", pad=3.0),
        )
    fig.savefig(
        "./curated-examples/figures/p-bit.png",
        dpi=300,
        bbox_inches="tight",
    )
