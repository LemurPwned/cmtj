"""
Smit-Beljers Static Equilibrium and FMR Frequency Calculation

GOAL:
This simulation uses the Smit-Beljers formalism to compute static equilibrium
magnetization orientations and ferromagnetic resonance (FMR) frequencies for
magnetic multilayers. Unlike dynamic LLG simulations, this approach finds
energy-minimized configurations and calculates resonance modes analytically,
providing fast characterization of magnetic properties and layer coupling.

CRITICAL SIMULATION MECHANISMS:
1. Energy Minimization: Finds equilibrium magnetization angles (theta, phi) by
   minimizing the total magnetic free energy including anisotropy, Zeeman,
   demagnetization, and exchange coupling terms
2. Hessian Analysis: Computes second derivatives of energy to determine stability
   and extract FMR frequencies from eigenvalues of the dynamic matrix
3. Spherical Coordinates: Uses (theta, phi) representation rather than Cartesian
   (mx, my, mz) for natural handling of energy landscapes on unit sphere
4. Field-Dependent Modes: Maps FMR frequencies as function of applied field to
   identify acoustic and optical modes in coupled multilayer systems
5. IEC Effects: Includes interlayer exchange coupling that splits degenerate
   modes and creates collective oscillation patterns
6. Multi-Layer Support: Handles arbitrary number of coupled layers with different
   magnetic properties, enabling SAF and spin-valve characterization

IMPORTANT UNITS NOTE (from AGENTS.md):
- Ms is in A/m for Smit-Beljers (NOT Tesla as in core Layer objects!)
- Ks (anisotropy) is in J/m³
- All fields are in A/m
- J (IEC) is in J/m²

The simulation generates field-swept FMR spectra showing resonance frequencies
and equilibrium angles. This enables rapid parameter extraction for complex
magnetic stacks without time-consuming full LLG integration.
"""

import contextlib

import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

from cmtj.models import LayerSB, Solver
from cmtj.utils import VectorObj, mu0

with contextlib.suppress(ImportError):
    import scienceplots  # noqa: F401

# CRITICAL: Smit-Beljers uses Ms in A/m, not T (different from core CMTJ!)
# Convert from Tesla to A/m: Ms[A/m] = Ms[T] / mu0
Ms1_T = 1.0  # T
Ms2_T = 0.9  # T
Ms1 = Ms1_T / mu0  # Convert to A/m for SB model
Ms2 = Ms2_T / mu0  # Convert to A/m for SB model

# Anisotropy constants in J/m³ (same units as core)
Ks1 = 400e3  # Perpendicular anisotropy layer 1
Ks2 = 350e3  # Perpendicular anisotropy layer 2

# In-plane anisotropy (Kv) - only phi angle matters, theta is ignored
# VectorObj(phi, theta, magnitude) for anisotropy axis
Kv_magnitude = 1e4  # Small in-plane anisotropy
phi_aniso = np.deg2rad(0)  # Along x-axis

# Create two coupled layers using LayerSB
layer1 = LayerSB(
    _id=0,
    thickness=1.5e-9,  # m
    Kv=VectorObj(phi_aniso, 0, Kv_magnitude),  # In-plane anisotropy
    Ks=Ks1,  # Perpendicular anisotropy
    Ms=Ms1,  # A/m (not T!)
    coordinate_system="spherical",  # SB uses spherical coordinates
)

layer2 = LayerSB(
    _id=1,
    thickness=1.2e-9,  # m
    Kv=VectorObj(phi_aniso, 0, Kv_magnitude),
    Ks=Ks2,
    Ms=Ms2,  # A/m (not T!)
    coordinate_system="spherical",
)

# IEC coupling between layers (J/m²)
# Following AGENTS.md: J range is ±0.001 to ±3.0 mJ/m²
J1 = -0.3e-3  # -0.3 mJ/m² (antiferromagnetic coupling)
J2 = 0.0  # No quadratic term

# Field sweep parameters
H_min = 0  # A/m
H_max = 500e3  # A/m (within ±0-±500e3 A/m range from AGENTS.md)
H_steps = 50

# Apply field in z direction (perpendicular to film)
theta_field = 0  # deg from z-axis
phi_field = 0  # deg in xy-plane

# Storage for results
field_values = np.linspace(H_min, H_max, H_steps)
frequencies = []
equilibrium_angles = []

print("Computing static equilibrium and FMR frequencies...")
print(f"Ms1 = {Ms1:.2e} A/m ({Ms1_T:.2f} T)")
print(f"Ms2 = {Ms2:.2e} A/m ({Ms2_T:.2f} T)")
print(f"J1 = {J1*1e3:.2f} mJ/m²")
print(f"Ks1 = {Ks1/1e3:.0f} kJ/m³")
print(f"Ks2 = {Ks2/1e3:.0f} kJ/m³\n")

# Initial guess for equilibrium position (theta1, phi1, theta2, phi2)
# Start with both layers pointing up (theta ≈ 0)
current_position = [0.1, 0.0, np.pi - 0.1, 0.0]  # Layer 1 up, Layer 2 down (AF coupling)

for H_magnitude in tqdm(field_values):
    # Create external field vector
    H_external = VectorObj.from_spherical(
        theta=np.deg2rad(theta_field),
        phi=np.deg2rad(phi_field),
        mag=H_magnitude,
    )

    # Create solver with both layers
    # J1 and J2 must be lists for coupling between layers
    solver = Solver(
        layers=[layer1, layer2],
        J1=[J1],  # List with coupling between layer1 and layer2
        J2=[J2],  # List with quadratic coupling
        H=H_external,
    )

    try:
        # Find equilibrium configuration and FMR frequencies
        # solve() returns (equilibrium_angles, frequencies, _optional_gradient_info)
        result = solver.solve(init_position=current_position)
        
        # Handle different return formats
        if len(result) == 3:
            eq_angles, freqs_ghz, _ = result
        else:
            eq_angles, freqs_ghz = result

        # Extract equilibrium angles (in radians)
        theta1, phi1, theta2, phi2 = eq_angles

        equilibrium_angles.append([theta1, phi1, theta2, phi2])

        # Frequencies are already in GHz from solve()
        frequencies.append(freqs_ghz)

        # Use current solution as initial guess for next iteration
        current_position = [theta1, phi1, theta2, phi2]

    except Exception as e:
        print(f"Warning: Failed at H={H_magnitude:.0f} A/m: {e}")
        equilibrium_angles.append([np.nan, np.nan, np.nan, np.nan])
        frequencies.append([np.nan, np.nan])

# Convert to arrays for plotting
equilibrium_angles = np.array(equilibrium_angles)
frequencies = np.array(frequencies)

# Create comprehensive plot
with plt.style.context(["science", "no-latex"]):
    fig, axes = plt.subplots(2, 2, figsize=(12, 10), dpi=300)

    # Plot 1: FMR frequencies vs field
    ax1 = axes[0, 0]
    H_ka = field_values / 1e3  # Convert to kA/m
    for i in range(frequencies.shape[1]):
        ax1.plot(
            H_ka,
            frequencies[:, i],
            linewidth=2,
            label=f"Mode {i+1}",
            marker="o",
            markersize=3,
        )
    ax1.set_xlabel("Applied Field (kA/m)")
    ax1.set_ylabel("FMR Frequency (GHz)")
    ax1.set_title("Ferromagnetic Resonance Spectrum")
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(bottom=0)

    # Plot 2: Equilibrium polar angles
    ax2 = axes[0, 1]
    ax2.plot(
        H_ka,
        np.rad2deg(equilibrium_angles[:, 0]),
        linewidth=2,
        color="crimson",
        label="Layer 1 θ",
    )
    ax2.plot(
        H_ka,
        np.rad2deg(equilibrium_angles[:, 2]),
        linewidth=2,
        color="navy",
        label="Layer 2 θ",
    )
    ax2.set_xlabel("Applied Field (kA/m)")
    ax2.set_ylabel("Polar Angle θ (degrees)")
    ax2.set_title("Equilibrium Magnetization (out-of-plane)")
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)

    # Plot 3: Equilibrium azimuthal angles
    ax3 = axes[1, 0]
    ax3.plot(
        H_ka,
        np.rad2deg(equilibrium_angles[:, 1]),
        linewidth=2,
        color="crimson",
        label="Layer 1 φ",
    )
    ax3.plot(
        H_ka,
        np.rad2deg(equilibrium_angles[:, 3]),
        linewidth=2,
        color="navy",
        label="Layer 2 φ",
    )
    ax3.set_xlabel("Applied Field (kA/m)")
    ax3.set_ylabel("Azimuthal Angle φ (degrees)")
    ax3.set_title("Equilibrium Magnetization (in-plane)")
    ax3.legend(fontsize=8)
    ax3.grid(True, alpha=0.3)

    # Plot 4: Mode splitting
    ax4 = axes[1, 1]
    if frequencies.shape[1] >= 2:
        mode_splitting = frequencies[:, 1] - frequencies[:, 0]
        ax4.plot(H_ka, mode_splitting, linewidth=2, color="forestgreen")
        ax4.set_xlabel("Applied Field (kA/m)")
        ax4.set_ylabel("Frequency Splitting (GHz)")
        ax4.set_title("Acoustic-Optical Mode Splitting")
        ax4.grid(True, alpha=0.3)
        ax4.axhline(y=0, color="red", linestyle="--", alpha=0.5)

    # Add parameter information
    param_text = (
        f"Parameters:\n"
        f"$M_{{s,1}}$ = {Ms1_T:.1f} T\n"
        f"$M_{{s,2}}$ = {Ms2_T:.1f} T\n"
        f"$K_{{s,1}}$ = {Ks1/1e3:.0f} kJ/m³\n"
        f"$K_{{s,2}}$ = {Ks2/1e3:.0f} kJ/m³\n"
        f"$J_1$ = {J1*1e3:.2f} mJ/m²\n"
        f"(Units: Ms in A/m for SB)"
    )
    fig.text(
        0.98,
        0.02,
        param_text,
        transform=fig.transFigure,
        fontsize=8,
        verticalalignment="bottom",
        horizontalalignment="right",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
    )

    fig.suptitle(
        "Smit-Beljers Model: Static Equilibrium & FMR", fontsize=14, y=0.995
    )
    fig.tight_layout(rect=[0, 0, 1, 0.99])
    fig.savefig(
        "./curated-examples/figures/sb-equilibrium.png",
        dpi=300,
        bbox_inches="tight",
    )

print(f"\nField sweep completed: {H_min/1e3:.0f} - {H_max/1e3:.0f} kA/m")
print(f"Number of modes found: {frequencies.shape[1]}")
print(f"FMR frequency at max field:")
for i in range(frequencies.shape[1]):
    print(f"  Mode {i+1}: {frequencies[-1, i]:.2f} GHz")
