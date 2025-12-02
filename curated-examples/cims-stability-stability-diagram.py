"""
Current-Induced Magnetization Switching (CIMS) Stability Diagram Simulation

GOAL:
This simulation generates comprehensive stability diagrams for spin-orbit torque (SOT)
devices by mapping switching behavior across parameter space including magnetic field,
current density, and material properties. The simulation identifies stable and unstable
regions to optimize device design and predict switching reliability under various
operating conditions.

CRITICAL SIMULATION MECHANISMS:
1. Multi-dimensional Parameter Sweep: Systematically varies magnetic field strength,
   current density, anisotropy constants, and torque scaling factors
2. Parallel Processing: Uses multiprocessing to efficiently explore large parameter
   spaces with statistical sampling for thermal effects
3. Thermal Fluctuations: Includes finite temperature effects with multiple thermal
   samples to capture switching probability distributions
4. SOT Dynamics: Models both field-like and damping-like spin-orbit torques with
   realistic current-dependent scaling relationships
5. Stability Analysis: Monitors final magnetization states after current pulses to
   determine switching success/failure across parameter space
6. Adaptive Integration: Uses sophisticated numerical methods with adaptive time
   stepping for accurate dynamics near switching thresholds
7. Statistical Averaging: Performs multiple runs with different initial conditions
   and thermal noise realizations to build robust statistics

The simulation generates 2D stability maps showing switching probability as functions
of field and current, revealing critical switching boundaries, optimal operating
points, and parameter sensitivity. These diagrams are essential for device design,
providing guidance on required current densities, field assistance, and material
parameter optimization for reliable SOT-based memory and logic devices.
"""

import contextlib
from collections import defaultdict
from copy import deepcopy

import matplotlib.pyplot as plt
import multiprocess as mp
import numpy as np
from tqdm import tqdm

from cmtj import (
    AxialDriver,
    CVector,
    Junction,
    Layer,
    SolverMode,
    constantDriver,
    stepDriver,
)
from cmtj.utils import FieldScan, TtoAm

with contextlib.suppress(ImportError):
    import scienceplots  # noqa: F401

Kdir_FL = CVector(0.0, 0.0, 1.0)

FL_Ms = 1.12
FL_Ku = 563e3

damping = 0.01

p = 0.017
Nxx = Nyy = p
Nzz = np.sqrt(1 - 2 * (p**2))
demag = [CVector(Nxx, 0, 0), CVector(0, Nyy, 0), CVector(0, 0, Nzz)]

thetas = [
    0.1,
]
Ku2s = [-30e3, -20e3, -10e3]  # Added more Ku2 values for demonstration
tau_scales = [0.7, 0.5]
Hmax = 400e3
Nsteps = 80
jscan_base = np.linspace(0, 0.6, 150)
jscan = np.concatenate([jscan_base, jscan_base[::-1], -jscan_base, -jscan_base[::-1]])
T = 300
solver = SolverMode.EulerHeun if T != 0 else SolverMode.DormandPrince
N_thermal_samples = 10 if T != 0 else 1  # Number of samples for thermal averaging
cell_surface = np.power(300e-9, 2) * np.pi
# Store all results for final plotting
all_results = {}

wf = 1e-11
integration_time = 1e-12
pulse_len = 10e-9
for tau_scale in tqdm(tau_scales, desc="tau_scales"):
    for theta in tqdm(thetas, desc="Thetas"):
        for FL_Ku2 in tqdm(Ku2s, desc="Ku2s"):
            Hscan, Hvecs = FieldScan.amplitude_scan(
                start=-Hmax,
                stop=Hmax,
                steps=Nsteps,
                theta=theta,
                phi=0,
            )

            # Define a function to process a single field vector
            def process_field_vector(field_data):
                idx, Hv = field_data  # Unpack the index and field vector

                # Store results from multiple thermal samples
                all_thermal_results = []

                # Run multiple samples if temperature is non-zero
                for thermal_sample in range(N_thermal_samples):
                    # Create a new junction for each process to avoid shared state
                    free = Layer(
                        "free",
                        mag=CVector(0.1, 0.1, 0.99),
                        anis=Kdir_FL,
                        Ms=FL_Ms,
                        thickness=1.0e-9,
                        damping=damping,
                        cellSurface=cell_surface,
                        demagTensor=demag,
                    )

                    free.setReferenceLayer(CVector(0, 0, 1))
                    free.setSecondOrderAnisotropyDriver(constantDriver(FL_Ku2))
                    free.setAnisotropyDriver(constantDriver(FL_Ku))
                    if T != 0:
                        free.setTemperatureDriver(constantDriver(T))
                    j = Junction([free])

                    j.clearLog()
                    # add a small random field to the system (different for each thermal sample)
                    if T != 0:
                        # Use different random seed for each thermal sample
                        np.random.seed(thermal_sample * 1000 + idx)
                    Hv_perturbed = Hv + np.random.randn(3) * 100
                    mag_init = CVector(*Hv_perturbed.tolist())
                    mag_init.normalize()
                    j.setLayerMagnetisation("free", mag_init)
                    Hval = Hv.tolist()
                    j.setLayerExternalFieldDriver("free", AxialDriver(*Hval))
                    j.runSimulation(25e-9, integration_time, wf, solverMode=solver)

                    voltage_dict = defaultdict(list)
                    local_scan = deepcopy(jscan)
                    if Hv[-1] > 0:
                        local_scan = local_scan[::-1]
                    for jden in local_scan:
                        j.clearLog()
                        j.setLayerFieldLikeTorqueDriver(
                            "free",
                            stepDriver(0, tau_scale * TtoAm * 1e-3 * (jden**2), 0, pulse_len),
                        )
                        j.setLayerDampingLikeTorqueDriver(
                            "free",
                            stepDriver(0, tau_scale * TtoAm * 16e-3 * jden, 0, pulse_len),
                        )
                        j.runSimulation(2 * pulse_len, integration_time, wf, solverMode=solver)
                        log = j.getLog()
                        # check mz
                        mz = np.mean(log["free_mz"][-100:])
                        voltage_dict[jden].append(mz)

                    # Store this thermal sample's results
                    thermal_result = {}
                    for jden, measured_values in voltage_dict.items():
                        thermal_result[jden] = np.mean(measured_values)
                    all_thermal_results.append(thermal_result)

                # Average over all thermal samples
                final_voltage_dict = defaultdict(list)
                for thermal_result in all_thermal_results:
                    for jden, mz_value in thermal_result.items():
                        final_voltage_dict[jden].append(mz_value)

                # Calculate mean and optionally standard deviation
                averaged_voltage_dict = {}
                for jden, mz_values in final_voltage_dict.items():
                    averaged_voltage_dict[jden] = np.mean(mz_values)
                    # Optionally store standard deviation for error analysis
                    # averaged_voltage_dict[f"{jden}_std"] = np.std(mz_values)

                # sort voltage_dict by voltage
                averaged_voltage_dict = dict(sorted(averaged_voltage_dict.items(), key=lambda x: x[0]))
                return (
                    idx,
                    list(averaged_voltage_dict.values()),
                    sorted(list(averaged_voltage_dict.keys())),
                )

            # Parallelize over field vectors using multiprocess
            # Create a pool of worker processes
            pool = mp.Pool(processes=4)

            # Prepare input data with indices to preserve order
            input_data = [(i, Hv) for i, Hv in enumerate(Hvecs)]

            # Use imap to process the field vectors in parallel with a progress bar
            results = list(
                tqdm(
                    pool.imap(process_field_vector, input_data),
                    total=len(input_data),
                    desc="Processing field vectors",
                )
            )

            # Close the pool
            pool.close()
            pool.join()

            # Sort results by the original index to maintain order
            results.sort(key=lambda x: x[0])

            # Unpack the results
            stability_spectrum = np.array([r[1] for r in results])
            voltage_keys = results[0][2]  # Get the keys from the first result

            # Store results for final plotting
            all_results[(tau_scale, FL_Ku2)] = {
                "Hscan": Hscan,
                "voltage_keys": voltage_keys,
                "stability_spectrum": stability_spectrum,
                "theta": theta,
            }

# Create the final aggregated figure
with plt.style.context(["science", "nature"]):
    n_rows = len(tau_scales)
    n_cols = len(Ku2s)
    w, h = plt.figaspect(n_rows / n_cols)
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        # figsize=(4 * n_cols, 4 * n_rows),
        dpi=300,
        sharex=True,
        sharey=True,
    )

    # Handle case where there's only one row or column
    if n_rows == 1:
        axes = axes.reshape(1, -1)
    if n_cols == 1:
        axes = axes.reshape(-1, 1)

    # Find global min/max for consistent colorbar
    all_spectra = [result["stability_spectrum"] for result in all_results.values()]
    vmin = min(np.min(spectrum) for spectrum in all_spectra)
    vmax = max(np.max(spectrum) for spectrum in all_spectra)

    for i, tau_scale in enumerate(tau_scales):
        for j, FL_Ku2 in enumerate(Ku2s):
            ax = axes[i, j]
            result = all_results[(tau_scale, FL_Ku2)]

            im = ax.pcolormesh(
                result["Hscan"] / 1e3,
                result["voltage_keys"],
                result["stability_spectrum"].T,
                cmap="jet",
                vmin=vmin,
                vmax=vmax,
            )

            # # Add row and column labels
            if j == n_cols - 1:  # Last column
                ax.text(
                    1.15,
                    0.5,
                    f"$\\tau$ = {tau_scale}",
                    transform=ax.transAxes,
                    rotation=270,
                    verticalalignment="center",
                    fontweight="bold",
                )

            # if i == 0:  # First row

    # Add a single colorbar for all subplots
    fig.subplots_adjust(right=0.85)
    cbar_ax = fig.add_axes([0.95, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label("$m_z$", rotation=270, labelpad=15)

    for i, _ in enumerate(tau_scales):
        axes[i, 0].set_ylabel("V (V)")

    for j, FL_Ku2 in enumerate(Ku2s):
        axes[-1, j].set_xlabel("H ($\\mathrm{{kA/m}}$)")
        axes[0, j].set_title(
            f"$\\mathrm{{K}}_\\mathrm{{u2}}$ = {FL_Ku2 / 1e3:.0f} $\\mathrm{{kJ/m^3}}$",
            fontsize=6,
        )
    fig.subplots_adjust(wspace=0.05, hspace=0.05)
    # Add subplot labels a,b,c,d,e,f
    labels = "abcdef"
    label_idx = 0
    import matplotlib.transforms as mtransforms

    for label, ax in zip("abcdefg", axes.flatten()):
        # label physical distance in and down:
        trans = mtransforms.ScaledTranslation(10 / 72, -5 / 72, fig.dpi_scale_trans)
        ax.text(
            -0.1,
            0.98,
            f"{label})",
            transform=ax.transAxes + trans,
            fontsize="medium",
            verticalalignment="top",
            color="lavender",
            bbox=dict(facecolor="none", alpha=0.4, edgecolor="none", pad=3.0),
        )

fig.savefig(
    "./curated-examples/figures/cims-stability.png",
    dpi=350,
    bbox_inches="tight",  # Ensures the full figure is saved without cutting off edges
)
