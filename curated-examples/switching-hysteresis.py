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

from cmtj import AxialDriver, CVector, Junction, Layer, ScalarDriver, constantDriver

with contextlib.suppress(ImportError):
    import scienceplots  # noqa: F401

# Layer parameters following AGENTS.md guidelines
# Ms in Tesla for core Layer objects
Ms = 1.0  # Saturation magnetization [T]

# Perpendicular anisotropy in typical range
Ku = 500e3  # PMA constant [J/m^3] - typical for PMA materials

# Anisotropy direction (perpendicular to film)
Kdir = CVector(0, 0, 1)

# Damping in recommended range (0.01-0.03)
damping = 0.02

# Demagnetization tensor (thin film - strong in-plane preference)
demag = [CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 1.0)]

# Create single magnetic layer
layer = Layer(
    "free",
    mag=CVector(0, 0, 1.0),  # Start with magnetization up
    anis=Kdir,
    Ms=Ms,
    thickness=1.5e-9,
    damping=damping,
    demagTensor=demag,
    cellSurface=np.pi * (40e-9) ** 2,  # Circular cross-section
)

# Set constant anisotropy
layer.setAnisotropyDriver(constantDriver(Ku))

# Create junction
junction = Junction([layer])

# Field sweep parameters (within ±0-±500e3 A/m range from AGENTS.md)
H_max = 300e3  # Maximum field [A/m]
H_min = -300e3  # Minimum field [A/m]
H_steps = 60  # Number of field steps

# Create field sweep: down sweep first, then up sweep
field_down = np.linspace(H_max, H_min, H_steps)
field_up = np.linspace(H_min, H_max, H_steps)
field_sweep = np.concatenate([field_down, field_up])

# Time parameters
# Following AGENTS.md: use dt=1e-12 for standard simulations
dt = 1e-12
relax_time = 5e-9  # 5 ns relaxation at each field step (within 1-500ns range)

# Storage for results
mz_values = []
applied_fields = []

print("Running hysteresis loop simulation...")
print(f"Ms = {Ms:.2f} T")
print(f"Ku = {Ku/1e3:.0f} kJ/m³")
print(f"Field range: {H_min/1e3:.0f} to {H_max/1e3:.0f} kA/m")
print(f"Damping α = {damping:.3f}\n")

# Run field sweep
for i, H_field in enumerate(field_sweep):
    # Set external field in z direction
    junction.setLayerExternalFieldDriver(
        "free",
        AxialDriver(
            ScalarDriver.getConstantDriver(0),
            ScalarDriver.getConstantDriver(0),
            ScalarDriver.getConstantDriver(H_field),
        ),
    )

    # Relax magnetization at this field
    junction.clearLog()
    junction.runSimulation(relax_time, dt, dt)

    # Get final magnetization state
    log = junction.getLog()
    mz_final = log["free_mz"][-1]

    mz_values.append(mz_final)
    applied_fields.append(H_field)

    # Progress indicator
    if (i + 1) % 20 == 0:
        print(f"Completed {i+1}/{len(field_sweep)} field points")

# Convert to arrays
mz_values = np.array(mz_values)
applied_fields = np.array(applied_fields)

# Split into down and up sweeps
split_idx = H_steps
H_down = applied_fields[:split_idx]
H_up = applied_fields[split_idx:]
mz_down = mz_values[:split_idx]
mz_up = mz_values[split_idx:]

# Find coercive fields (where mz crosses zero)
# Down sweep: find where mz goes from positive to negative
try:
    idx_down = np.where(np.diff(np.sign(mz_down)))[0][0]
    Hc_down = H_down[idx_down]
except IndexError:
    Hc_down = np.nan

# Up sweep: find where mz goes from negative to positive
try:
    idx_up = np.where(np.diff(np.sign(mz_up)))[0][0]
    Hc_up = H_up[idx_up]
except IndexError:
    Hc_up = np.nan

print(f"\nCoercive field (down sweep): {Hc_down/1e3:.1f} kA/m")
print(f"Coercive field (up sweep): {Hc_up/1e3:.1f} kA/m")
print(f"Coercivity asymmetry: {abs(Hc_down - Hc_up)/1e3:.1f} kA/m")

# Create simplified M(H) plot
with plt.style.context(["science", "no-latex"]):
    fig, ax = plt.subplots(figsize=(6, 5), dpi=300)

    # Main hysteresis loop (mz vs H)
    ax.plot(H_down / 1e3, mz_down, "-", color="crimson", linewidth=2, label="Down sweep")
    ax.plot(H_up / 1e3, mz_up, "-", color="navy", linewidth=2, label="Up sweep")
    ax.axhline(y=0, color="gray", linestyle="--", alpha=0.5)
    ax.axvline(x=0, color="gray", linestyle="--", alpha=0.5)
    if not np.isnan(Hc_down):
        ax.axvline(x=Hc_down / 1e3, color="red", linestyle=":", alpha=0.7, label=f"$H_c$↓={Hc_down/1e3:.1f} kA/m")
    if not np.isnan(Hc_up):
        ax.axvline(x=Hc_up / 1e3, color="blue", linestyle=":", alpha=0.7, label=f"$H_c$↑={Hc_up/1e3:.1f} kA/m")
    ax.set_xlabel(r"$H$ (kA/m)")
    ax.set_ylabel(r"$m_z$")
    ax.set_title(f"$M_s$ = {Ms:.1f} T, $K_u$ = {Ku/1e3:.0f} kJ/m³, $\\alpha$ = {damping:.3f}")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_ylim([-1.1, 1.1])

    fig.tight_layout()
    fig.savefig(
        "./curated-examples/figures/switching-hysteresis.png",
        dpi=300,
        bbox_inches="tight",
    )

print("\nSimulation complete!")
print(f"Hysteresis loop figure saved to curated-examples/figures/switching-hysteresis.png")
