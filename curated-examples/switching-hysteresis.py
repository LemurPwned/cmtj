"""
Magnetization Switching and Hysteresis Loop Simulation

GOAL:
This simulation demonstrates basic magnetization switching dynamics in a single
ferromagnetic layer under an applied magnetic field sweep. It captures the hysteretic
behavior characteristic of magnetic memory elements and illustrates the fundamental
relationship between applied field and magnetization state.

CRITICAL SIMULATION MECHANISMS:
1. Single Layer Dynamics: Models a free magnetic layer with perpendicular magnetic
   anisotropy (PMA) that can switch between up and down states
2. Field Sweep Protocol: Applies ramped external magnetic field to drive magnetization
   reversal through coherent rotation or domain nucleation/propagation
3. Hysteresis Behavior: Captures history-dependent switching with distinct coercive
   fields for up-to-down and down-to-up transitions
4. Dynamic Switching: Uses full LLG equation integration to capture realistic
   switching trajectories including precession and damping effects
5. Anisotropy Energy: Includes uniaxial perpendicular anisotropy that creates
   energy barriers and defines the easy axis for magnetization
6. Demagnetization: Accounts for shape anisotropy that favors in-plane orientation
   for thin films, competing with perpendicular anisotropy

The simulation generates magnetization vs field (M-H) hysteresis loops and
time-domain switching trajectories. This enables characterization of coercive
field, switching field distribution, and energy barrier heights - critical
parameters for magnetic memory (MRAM) and sensor applications.
"""

import contextlib

import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

from cmtj import AxialDriver, CVector, Junction, Layer, ScalarDriver, constantDriver

with contextlib.suppress(ImportError):
    import scienceplots  # noqa: F401

# Layer parameters from documentation (docs/experimental-methods/examples.ipynb)
# Ms in Tesla for core Layer objects
Ms = 1.03  # Saturation magnetization [T] - from examples.ipynb

# In-plane anisotropy. Use a stronger easy-axis term so the example shows a
# clear finite-field hysteresis loop instead of near-reversible rotation.
Ku = 8.0e3  # Anisotropy constant [J/m^3]

# Anisotropy direction - IN-PLANE (along x) for in-plane hysteresis
Kdir = CVector(1, 0, 0)

# Damping from documentation
damping = 0.024  # Gilbert damping - from examples.ipynb

# Standard thin film demagnetization tensor from documentation
demag = [CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 1.0)]

# Create single magnetic layer with in-plane initial magnetization
# Start slightly off-axis to allow switching
layer = Layer(
    "free",
    mag=CVector(1.0, 0.1, 0.1),  # Start near +x with small perturbation
    anis=Kdir,
    Ms=Ms,
    thickness=2.1e-9,  # From examples.ipynb
    damping=damping,
    demagTensor=demag,
    cellSurface=15e-9 * 15e-9 * np.pi,  # From examples.ipynb
)

# Set constant anisotropy
layer.setAnisotropyDriver(constantDriver(Ku))

# Create junction
junction = Junction([layer])

# Field sweep parameters - IN-PLANE field along easy axis (x) for switching
H_max = 70e3  # Maximum field [A/m] in x direction
H_min = -70e3  # Minimum field [A/m]
H_steps = 120  # Number of field steps
H_bias_y = 100.0  # tiny transverse bias to avoid getting stuck on the easy axis

# Create field sweep: 0 -> Hmin -> Hmax -> 0.
# This makes the negative-going and positive-going switching branches explicit.
field_to_negative = np.linspace(0, H_min, H_steps, endpoint=False)
field_negative_to_positive = np.linspace(H_min, H_max, 2 * H_steps, endpoint=False)
field_positive_to_zero = np.linspace(H_max, 0, H_steps)
field_sweep = np.concatenate([field_to_negative, field_negative_to_positive, field_positive_to_zero])

# Time parameters
dt = 1e-13
relax_time = 10e-9

# Storage for results - track mx since field is in x direction
mx_values = []
applied_fields = []

print("Running hysteresis loop simulation...")
print(f"Ms = {Ms:.2f} T")
print(f"Ku = {Ku / 1e3:.1f} kJ/m³")
print(f"Field range: {H_min / 1e3:.0f} to {H_max / 1e3:.0f} kA/m")
print(f"Damping α = {damping:.3f}\n")

# Prepare the layer in a reproducible positive state before recording the loop.
junction.setLayerExternalFieldDriver(
    "free",
    AxialDriver(
        ScalarDriver.getConstantDriver(H_max),
        ScalarDriver.getConstantDriver(H_bias_y),
        ScalarDriver.getConstantDriver(0),
    ),
)
junction.clearLog()
junction.runSimulation(relax_time, dt, dt)

# Run field sweep
for H_field in tqdm(field_sweep):
    # Set external field in X direction (along easy axis)
    junction.setLayerExternalFieldDriver(
        "free",
        AxialDriver(
            ScalarDriver.getConstantDriver(H_field),
            ScalarDriver.getConstantDriver(H_bias_y),
            ScalarDriver.getConstantDriver(0),
        ),
    )

    # Relax magnetization at this field
    junction.clearLog()
    junction.runSimulation(relax_time, dt, dt)

    # Get final magnetization state
    log = junction.getLog()
    mx_final = log["free_mx"][-1]

    mx_values.append(mx_final)
    applied_fields.append(H_field)

# Convert to arrays
mx_values = np.array(mx_values)
applied_fields = np.array(applied_fields)

# Split into the three sweep branches
first_split = len(field_to_negative)
second_split = first_split + len(field_negative_to_positive)

H_down = applied_fields[:first_split]
mx_down = mx_values[:first_split]

H_up = applied_fields[first_split:second_split]
mx_up = mx_values[first_split:second_split]

H_return = applied_fields[second_split:]
mx_return = mx_values[second_split:]

# Find coercive fields (where mx crosses zero)
# Down sweep: find where mx goes from positive to negative
try:
    idx_down = np.where(np.diff(np.sign(mx_down)))[0][0]
    Hc_down = H_down[idx_down]
except IndexError:
    Hc_down = np.nan

# Up sweep: find where mx goes from negative to positive
try:
    idx_up = np.where(np.diff(np.sign(mx_up)))[0][0]
    Hc_up = H_up[idx_up]
except IndexError:
    Hc_up = np.nan

print(f"\nCoercive field (down sweep): {Hc_down / 1e3:.1f} kA/m")
print(f"Coercive field (up sweep): {Hc_up / 1e3:.1f} kA/m")
print(f"Coercivity asymmetry: {abs(Hc_down - Hc_up) / 1e3:.1f} kA/m")

# Create simplified M(H) plot
with plt.style.context(["science", "no-latex"]):
    fig, ax = plt.subplots(figsize=(6, 5), dpi=300)

    # Main hysteresis loop (mx vs H)
    ax.plot(H_down / 1e3, mx_down, "-", color="crimson", linewidth=2, label="0 → $H_{min}$")
    ax.plot(H_up / 1e3, mx_up, "-", color="navy", linewidth=2, label="$H_{min}$ → $H_{max}$")
    ax.plot(H_return / 1e3, mx_return, "--", color="darkgreen", linewidth=1.8, label="$H_{max}$ → 0")
    ax.axhline(y=0, color="gray", linestyle="--", alpha=0.5)
    ax.axvline(x=0, color="gray", linestyle="--", alpha=0.5)
    if not np.isnan(Hc_down):
        ax.axvline(x=Hc_down / 1e3, color="red", linestyle=":", alpha=0.7, label=f"$H_c$↓={Hc_down / 1e3:.1f} kA/m")
    if not np.isnan(Hc_up):
        ax.axvline(x=Hc_up / 1e3, color="blue", linestyle=":", alpha=0.7, label=f"$H_c$↑={Hc_up / 1e3:.1f} kA/m")
    ax.set_xlabel(r"$H$ (kA/m)")
    ax.set_ylabel(r"$m_x$")
    ax.set_title(f"$M_s$ = {Ms:.2f} T, $K_u$ = {Ku / 1e3:.1f} kJ/m³, $\\alpha$ = {damping:.3f}")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_ylim([-1.1, 1.1])

    fig.tight_layout()
    fig.savefig(
        "./curated-examples/figures/switching-hysteresis.png",
        dpi=300,
        bbox_inches="tight",
    )
