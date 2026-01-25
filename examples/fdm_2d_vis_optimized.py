"""
FDM 2D Visualization Example (Optimized)
=========================================

This example demonstrates FDM simulation with 2D spatial plots of the
magnetization component mz over a 20x20x1 grid.

Optimized version that runs efficiently while still showing dynamics.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import cmtj
from cmtj.utils import get_full_demag_tensor, convert_tensor_to_cpp
import time

# Physical parameters
Ms = 1e6  # Saturation magnetization [A/m]
Ku = 1e4  # Anisotropy constant [J/m^3]
damping = 0.05  # Gilbert damping
thickness = 2e-9  # Layer thickness [m]

# Grid parameters (20x20x1)
width = 100e-9  # Total width in x [m]
length = 100e-9  # Total length in y [m]
cellSizeXY = 5e-9  # Cell size -> 20x20 grid
nz = 1

# Simulation parameters - using dt=1e-14 as requested
totalTime = 5e-9  # 5 nanoseconds total
timeStep = 1e-14  # 0.01 ps time step
n_segments = 10  # Break into segments to capture state
time_per_segment = totalTime / n_segments

print("=" * 70)
print("FDM 2D Visualization Simulation")
print("=" * 70)
print(f"Grid: 20 x 20 x 1")
print(f"Cell size: {cellSizeXY*1e9:.1f}nm")
print(f"Time step: {timeStep*1e12:.2f} ps (1e-14 s)")
print(f"Total time: {totalTime*1e9:.1f} ns")
print(f"Total time steps: {int(totalTime/timeStep):,}")
print(f"Segments: {n_segments} (capturing {n_segments+1} snapshots)")
print(f"Steps per segment: {int(time_per_segment/timeStep):,}")
print("=" * 70)

# Create layer template
layer_id = "free_layer"
initial_mag = cmtj.CVector(0, 0, 1)
anisotropy_axis = cmtj.CVector(0, 0, 1)

demagTensor = [
    cmtj.CVector(0, 0, 0),
    cmtj.CVector(0, 0, 0),
    cmtj.CVector(0, 0, 1)
]

layer_template = cmtj.Layer(
    id=layer_id,
    mag=initial_mag,
    anis=anisotropy_axis,
    Ms=Ms,
    thickness=thickness,
    cellSurface=cellSizeXY * cellSizeXY,
    demagTensor=demagTensor,
    damping=damping
)

layer_template.setAnisotropyDriver(cmtj.constantDriver(Ku))

# Create FDM grid
grid_spec = cmtj.fdm.FDMGridSpec(
    width=width,
    length=length,
    thickness=thickness,
    cellSizeXY=cellSizeXY,
    nz=nz
)

print(f"\nGrid: nx={grid_spec.nx}, ny={grid_spec.ny}, nz={grid_spec.nz}")
print(f"Cell dimensions: {grid_spec.dx*1e9:.2f}nm x {grid_spec.dy*1e9:.2f}nm x {grid_spec.dz*1e9:.2f}nm")

# Create FDM layer
fdm_layer = cmtj.fdm.FDMLayer(layer_template, grid_spec)

# Compute demag tensor
print("\nComputing demagnetization tensor...")
start_demag = time.time()
n = (grid_spec.nx, grid_spec.ny, grid_spec.nz)
dx_tuple = (grid_spec.dx, grid_spec.dy, grid_spec.dz)
demag_tensor = get_full_demag_tensor(n, dx_tuple)
cpp_demag_tensor = convert_tensor_to_cpp(demag_tensor)
fdm_layer.setDemagTensor(cpp_demag_tensor)
print(f"Demag tensor computed and set in {time.time()-start_demag:.2f}s")

# Set initial magnetization - vortex pattern
print("\nSetting vortex initial state...")
mags = []
center_x = grid_spec.nx / 2
center_y = grid_spec.ny / 2

for iz in range(grid_spec.nz):
    for iy in range(grid_spec.ny):
        for ix in range(grid_spec.nx):
            dx_cell = ix - center_x
            dy_cell = iy - center_y
            r = np.sqrt(dx_cell**2 + dy_cell**2)

            if r < 5:
                # Central vortex core
                angle = np.arctan2(dy_cell, dx_cell)
                mx = 0.4 * np.cos(angle)
                my = 0.4 * np.sin(angle)
                mz = 0.8
            else:
                # Outer region with small perturbation
                mx = 0.02 * (np.random.rand() - 0.5)
                my = 0.02 * (np.random.rand() - 0.5)
                mz = 1.0

            norm = np.sqrt(mx**2 + my**2 + mz**2)
            mags.append(cmtj.CVector(mx/norm, my/norm, mz/norm))

fdm_layer.setMagnetisationGrid(mags)

# Create junction
junction = cmtj.fdm.FDMJunction([fdm_layer])

# Apply external field
H_ext = 8e4  # External field [A/m]
field_driver = cmtj.AxialDriver(
    cmtj.constantDriver(H_ext),
    cmtj.NullDriver(),
    cmtj.constantDriver(H_ext * 0.2)  # Small z component
)
junction.setLayerExternalFieldDriver(layer_id, field_driver)

# Run simulation in segments
print(f"\nRunning simulation ({n_segments} segments)...")
print("This will take a few minutes...")

snapshots = []
snapshot_times = []

# Capture initial state
snapshots.append(junction.getLayerMagnetisationGrid(layer_id))
snapshot_times.append(0.0)

for i in range(n_segments):
    segment_start = time.time()

    junction.runSimulation(
        totalTime=time_per_segment,
        timeStep=timeStep,
        writeFrequency=time_per_segment,
        verbose=False,
        solverMode=cmtj.RK4
    )

    # Capture state
    snapshots.append(junction.getLayerMagnetisationGrid(layer_id))
    snapshot_times.append((i + 1) * time_per_segment)

    elapsed = time.time() - segment_start
    print(f"  Segment {i+1}/{n_segments} complete ({elapsed:.1f}s, t={snapshot_times[-1]*1e9:.2f}ns)")

print("\nSimulation complete!")

# Convert to numpy arrays
def mags_to_array(mags, nx, ny):
    arr = np.zeros((nx, ny, 3))
    idx = 0
    for iy in range(ny):
        for ix in range(nx):
            arr[ix, iy, 0] = mags[idx].x
            arr[ix, iy, 1] = mags[idx].y
            arr[ix, iy, 2] = mags[idx].z
            idx += 1
    return arr

# Create 2D plots
print("\nGenerating 2D plots of mz component...")

fig, axes = plt.subplots(2, 3, figsize=(16, 11))
axes = axes.flatten()

# Plot 6 snapshots
plot_indices = [0, 2, 4, 6, 8, 10]

for i, idx in enumerate(plot_indices):
    mag_array = mags_to_array(snapshots[idx], grid_spec.nx, grid_spec.ny)
    mz = mag_array[:, :, 2]
    mx = mag_array[:, :, 0]
    my = mag_array[:, :, 1]

    ax = axes[i]

    # Plot mz as heatmap
    im = ax.imshow(mz.T, cmap='RdBu_r', vmin=-1, vmax=1, origin='lower',
                   extent=[0, width*1e9, 0, length*1e9], interpolation='bilinear')

    # Overlay in-plane magnetization vectors
    step = 2  # Show every 2nd vector
    x_pos = np.arange(0, grid_spec.nx, step) * grid_spec.dx * 1e9 + grid_spec.dx * 1e9 / 2
    y_pos = np.arange(0, grid_spec.ny, step) * grid_spec.dy * 1e9 + grid_spec.dy * 1e9 / 2
    X, Y = np.meshgrid(x_pos, y_pos)

    U = mx[::step, ::step].T
    V = my[::step, ::step].T

    ax.quiver(X, Y, U, V, color='white', alpha=0.7, scale=12, width=0.004,
              edgecolor='black', linewidth=0.5)

    ax.set_title(f't = {snapshot_times[idx]*1e9:.2f} ns', fontsize=12, fontweight='bold')
    ax.set_xlabel('x (nm)', fontsize=10)
    ax.set_ylabel('y (nm)', fontsize=10)
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('$m_z$', fontsize=10)

plt.suptitle(f'FDM Simulation: $m_z$ Component Evolution\\n' +
             f'20×20×1 grid, dt={timeStep*1e12:.2f}ps, H_ext={H_ext/1e3:.0f}kA/m',
             fontsize=14, fontweight='bold')
plt.tight_layout()

# Save
output_file = 'fdm_2d_mz_snapshots.png'
plt.savefig(output_file, dpi=200, bbox_inches='tight')
print(f"Plot saved: {output_file}")

# Create simple animation
print("\nGenerating animation...")
fig_anim, ax_anim = plt.subplots(figsize=(10, 9))

def update_frame(frame_idx):
    ax_anim.clear()
    mag_array = mags_to_array(snapshots[frame_idx], grid_spec.nx, grid_spec.ny)
    mz = mag_array[:, :, 2]
    mx = mag_array[:, :, 0]
    my = mag_array[:, :, 1]

    im = ax_anim.imshow(mz.T, cmap='RdBu_r', vmin=-1, vmax=1, origin='lower',
                        extent=[0, width*1e9, 0, length*1e9], interpolation='bilinear')

    # Quiver
    step = 2
    x_pos = np.arange(0, grid_spec.nx, step) * grid_spec.dx * 1e9 + grid_spec.dx * 1e9 / 2
    y_pos = np.arange(0, grid_spec.ny, step) * grid_spec.dy * 1e9 + grid_spec.dy * 1e9 / 2
    X, Y = np.meshgrid(x_pos, y_pos)
    U = mx[::step, ::step].T
    V = my[::step, ::step].T

    ax_anim.quiver(X, Y, U, V, color='white', alpha=0.7, scale=12, width=0.004,
                   edgecolor='black', linewidth=0.5)

    ax_anim.set_title(f'FDM Magnetization: t = {snapshot_times[frame_idx]*1e9:.2f} ns',
                      fontsize=14, fontweight='bold')
    ax_anim.set_xlabel('x (nm)', fontsize=12)
    ax_anim.set_ylabel('y (nm)', fontsize=12)
    ax_anim.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
    return [im]

anim = FuncAnimation(fig_anim, update_frame, frames=len(snapshots),
                     interval=200, blit=False, repeat=True)

anim_file = 'fdm_2d_mz_animation.gif'
anim.save(anim_file, writer='pillow', fps=5, dpi=100)
print(f"Animation saved: {anim_file}")

print("\n" + "=" * 70)
print("Visualization complete!")
print(f"  - Static plots: {output_file}")
print(f"  - Animation: {anim_file}")
print("=" * 70)
