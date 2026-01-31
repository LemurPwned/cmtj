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
mx_values = []
my_values = []
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
    mx_final = log["free_mx"][-1]
    my_final = log["free_my"][-1]

    mz_values.append(mz_final)
    mx_values.append(mx_final)
    my_values.append(my_final)
    applied_fields.append(H_field)

    # Progress indicator
    if (i + 1) % 20 == 0:
        print(f"Completed {i+1}/{len(field_sweep)} field points")

# Convert to arrays
mz_values = np.array(mz_values)
mx_values = np.array(mx_values)
my_values = np.array(my_values)
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

# Create plots
with plt.style.context(["science", "no-latex"]):
    fig, axes = plt.subplots(2, 2, figsize=(12, 10), dpi=300)

    # Plot 1: Main hysteresis loop (mz vs H)
    ax1 = axes[0, 0]
    ax1.plot(H_down / 1e3, mz_down, "o-", color="crimson", linewidth=2, markersize=4, label="Down sweep")
    ax1.plot(H_up / 1e3, mz_up, "s-", color="navy", linewidth=2, markersize=4, label="Up sweep")
    ax1.axhline(y=0, color="gray", linestyle="--", alpha=0.5)
    ax1.axvline(x=0, color="gray", linestyle="--", alpha=0.5)
    if not np.isnan(Hc_down):
        ax1.axvline(x=Hc_down / 1e3, color="red", linestyle=":", alpha=0.7, label=f"Hc↓={Hc_down/1e3:.1f}")
    if not np.isnan(Hc_up):
        ax1.axvline(x=Hc_up / 1e3, color="blue", linestyle=":", alpha=0.7, label=f"Hc↑={Hc_up/1e3:.1f}")
    ax1.set_xlabel("Applied Field (kA/m)")
    ax1.set_ylabel("$m_z$ (normalized)")
    ax1.set_title("Hysteresis Loop (Out-of-Plane)")
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim([-1.1, 1.1])

    # Plot 2: In-plane components
    ax2 = axes[0, 1]
    ax2.plot(applied_fields / 1e3, mx_values, color="forestgreen", linewidth=2, label="$m_x$")
    ax2.plot(applied_fields / 1e3, my_values, color="orange", linewidth=2, label="$m_y$")
    ax2.set_xlabel("Applied Field (kA/m)")
    ax2.set_ylabel("In-plane magnetization")
    ax2.set_title("In-Plane Components")
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)

    # Plot 3: Phase space trajectory
    ax3 = axes[1, 0]
    # Color by field value
    colors = applied_fields / 1e3
    scatter = ax3.scatter(mx_values, my_values, c=colors, s=30, cmap="viridis", alpha=0.7)
    ax3.scatter([mx_values[0]], [my_values[0]], color="green", s=100, marker="o", zorder=10, label="Start")
    ax3.scatter([mx_values[-1]], [my_values[-1]], color="red", s=100, marker="s", zorder=10, label="End")
    ax3.set_xlabel("$m_x$")
    ax3.set_ylabel("$m_y$")
    ax3.set_title("In-Plane Phase Space")
    ax3.legend(fontsize=8)
    ax3.grid(True, alpha=0.3)
    ax3.set_aspect("equal")
    cbar = plt.colorbar(scatter, ax=ax3)
    cbar.set_label("H (kA/m)", fontsize=8)

    # Plot 4: Total magnetization magnitude
    ax4 = axes[1, 1]
    m_total = np.sqrt(mx_values**2 + my_values**2 + mz_values**2)
    ax4.plot(applied_fields / 1e3, m_total, color="purple", linewidth=2)
    ax4.axhline(y=1.0, color="red", linestyle="--", alpha=0.5, label="Ideal (|m|=1)")
    ax4.set_xlabel("Applied Field (kA/m)")
    ax4.set_ylabel("Total |m|")
    ax4.set_title("Magnetization Magnitude Conservation")
    ax4.legend(fontsize=8)
    ax4.grid(True, alpha=0.3)
    ax4.set_ylim([0.995, 1.005])

    # Add parameter information
    param_text = (
        f"Parameters:\n"
        f"$M_s$ = {Ms:.1f} T\n"
        f"$K_u$ = {Ku/1e3:.0f} kJ/m³\n"
        f"$\\alpha$ = {damping:.3f}\n"
        f"Thickness = {layer.thickness*1e9:.1f} nm\n"
        f"Relax time = {relax_time*1e9:.1f} ns"
    )
    fig.text(
        0.98,
        0.02,
        param_text,
        transform=fig.transFigure,
        fontsize=8,
        verticalalignment="bottom",
        horizontalalignment="right",
        bbox=dict(boxstyle="round", facecolor="lightgreen", alpha=0.5),
    )

    fig.suptitle("Magnetization Switching and Hysteresis", fontsize=14, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.99])
    fig.savefig(
        "./curated-examples/figures/switching-hysteresis.png",
        dpi=300,
        bbox_inches="tight",
    )

print("\nSimulation complete!")
print(f"Hysteresis loop figure saved to curated-examples/figures/switching-hysteresis.png")
