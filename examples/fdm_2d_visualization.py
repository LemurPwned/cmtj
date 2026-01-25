"""
FDM 2D Visualization Example
=============================

This example demonstrates FDM simulation with 2D spatial plots of the
magnetization component mz over a 20x20x1 grid.

We simulate magnetization dynamics under an external field and visualize
the spatial distribution of mz at different time snapshots.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import cmtj
from cmtj.utils import get_full_demag_tensor, convert_tensor_to_cpp

# Physical parameters
Ms = 1e6  # Saturation magnetization [A/m]
Ku = 1e4  # Anisotropy constant [J/m^3]
damping = 0.05  # Gilbert damping (increased for faster relaxation)
thickness = 2e-9  # Layer thickness [m]

# Grid parameters (20x20x1)
width = 100e-9  # Total width in x [m]
length = 100e-9  # Total length in y [m]
cellSizeXY = 5e-9  # Cell size in xy plane [m] -> 20x20 grid
nz = 1  # Single layer

# Simulation parameters
totalTime = 1e-9  # 1 nanosecond (reduced for reasonable runtime with dt=1e-14)
timeStep = 1e-14  # 0.01 picoseconds (as requested)
writeFrequency = 5e-11  # Write every 0.05 ns -> 20 snapshots

print("=" * 60)
print("FDM 2D Visualization Simulation")
print("=" * 60)
print(f"Grid: 20 x 20 x 1")
print(f"Cell size: {cellSizeXY*1e9:.1f}nm")
print(f"Time step: {timeStep*1e12:.2f} ps")
print(f"Total time: {totalTime*1e9:.1f} ns")
print(f"Write frequency: {writeFrequency*1e9:.2f} ns")
print(f"Expected snapshots: {int(totalTime/writeFrequency)}")
print("=" * 60)

# Create a base Layer template
layer_id = "free_layer"
initial_mag = cmtj.CVector(0, 0, 1)  # Initial magnetization along +z
anisotropy_axis = cmtj.CVector(0, 0, 1)  # Anisotropy along +z

# DemagTensor for Layer (placeholder, will use FDM per-cell tensors)
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

# Set anisotropy
layer_template.setAnisotropyDriver(cmtj.constantDriver(Ku))

# Create FDM grid specification
grid_spec = cmtj.fdm.FDMGridSpec(
    width=width,
    length=length,
    thickness=thickness,
    cellSizeXY=cellSizeXY,
    nz=nz
)

print(f"\nGrid verification:")
print(f"  nx={grid_spec.nx}, ny={grid_spec.ny}, nz={grid_spec.nz}")
print(f"  dx={grid_spec.dx*1e9:.2f}nm, dy={grid_spec.dy*1e9:.2f}nm, dz={grid_spec.dz*1e9:.2f}nm")
print(f"  Total cells: {grid_spec.nx * grid_spec.ny * grid_spec.nz}")

# Create FDM layer
fdm_layer = cmtj.fdm.FDMLayer(layer_template, grid_spec)

# Compute and set demagnetization tensor
print("\nComputing demagnetization tensor...")
n = (grid_spec.nx, grid_spec.ny, grid_spec.nz)
dx = (grid_spec.dx, grid_spec.dy, grid_spec.dz)
demag_tensor = get_full_demag_tensor(n, dx)
print(f"Demag tensor shape: {demag_tensor.shape}")

cpp_demag_tensor = convert_tensor_to_cpp(demag_tensor)
fdm_layer.setDemagTensor(cpp_demag_tensor)
print("Demagnetization tensor set successfully")

# Set initial magnetization with vortex-like pattern
print("\nSetting initial vortex-like magnetization pattern...")
mags = []
center_x = grid_spec.nx / 2
center_y = grid_spec.ny / 2

for iz in range(grid_spec.nz):
    for iy in range(grid_spec.ny):
        for ix in range(grid_spec.nx):
            # Distance from center
            dx_cell = ix - center_x
            dy_cell = iy - center_y
            r = np.sqrt(dx_cell**2 + dy_cell**2)

            # Create a vortex-like initial state
            # Outer region: pointing up (+z)
            # Center region: tilted based on angle
            if r < 5:  # Central 5-cell radius
                angle = np.arctan2(dy_cell, dx_cell)
                # Tilt magnetization in-plane at center
                mx = 0.3 * np.cos(angle)
                my = 0.3 * np.sin(angle)
                mz = 0.9
            else:
                # Outer region: mostly +z with small random perturbation
                mx = 0.05 * (np.random.rand() - 0.5)
                my = 0.05 * (np.random.rand() - 0.5)
                mz = 1.0

            # Normalize
            norm = np.sqrt(mx**2 + my**2 + mz**2)
            mags.append(cmtj.CVector(mx/norm, my/norm, mz/norm))

fdm_layer.setMagnetisationGrid(mags)
print(f"Set initial magnetization for {len(mags)} cells")

# Create FDM junction
junction = cmtj.fdm.FDMJunction([fdm_layer])

# Apply external field in +x direction to drive dynamics
H_ext = 5e4  # Moderate external field [A/m]
field_driver = cmtj.AxialDriver(
    cmtj.constantDriver(H_ext),  # Field in x direction
    cmtj.NullDriver(),
    cmtj.NullDriver()
)
junction.setLayerExternalFieldDriver(layer_id, field_driver)

print("\nRunning FDM simulation...")
print(f"This may take a while (500,000 time steps)...")

# Store magnetization snapshots during simulation
# We'll need to modify the approach since runSimulation doesn't return intermediate states
# For now, we'll run the simulation and then create plots

# Since we can't get intermediate snapshots easily, let's run shorter simulations
# and capture the state at each interval
snapshots = []
snapshot_times = []
n_snapshots = 20  # Reduced for faster execution
time_per_snapshot = totalTime / n_snapshots

print(f"\nRunning simulation in {n_snapshots} segments to capture snapshots...")

for i in range(n_snapshots + 1):
    # Get current magnetization
    current_mags = junction.getLayerMagnetisationGrid(layer_id)
    snapshots.append(current_mags)
    snapshot_times.append(i * time_per_snapshot)

    if i < n_snapshots:
        # Run for one time segment
        junction.runSimulation(
            totalTime=time_per_snapshot,
            timeStep=timeStep,
            writeFrequency=time_per_snapshot,  # Write at end only
            verbose=False,
            solverMode=cmtj.RK4
        )
        if (i + 1) % 10 == 0:
            print(f"  Progress: {i+1}/{n_snapshots} segments completed")

print("\nSimulation complete!")

# Convert snapshots to numpy arrays for plotting
def mags_to_array(mags, nx, ny):
    """Convert list of CVector magnetizations to numpy array (nx, ny, 3)"""
    arr = np.zeros((nx, ny, 3))
    idx = 0
    for iy in range(ny):
        for ix in range(nx):
            arr[ix, iy, 0] = mags[idx].x
            arr[ix, iy, 1] = mags[idx].y
            arr[ix, iy, 2] = mags[idx].z
            idx += 1
    return arr

# Create visualization
print("\nGenerating 2D plots...")

# Select 6 time points to plot
plot_indices = [0, 4, 8, 12, 16, 20]
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
axes = axes.flatten()

for i, idx in enumerate(plot_indices):
    mag_array = mags_to_array(snapshots[idx], grid_spec.nx, grid_spec.ny)
    mz = mag_array[:, :, 2]

    ax = axes[i]
    im = ax.imshow(mz.T, cmap='RdBu_r', vmin=-1, vmax=1, origin='lower',
                   extent=[0, width*1e9, 0, length*1e9])
    ax.set_title(f't = {snapshot_times[idx]*1e9:.2f} ns')
    ax.set_xlabel('x (nm)')
    ax.set_ylabel('y (nm)')

    # Add quiver plot to show in-plane components
    mx = mag_array[:, :, 0]
    my = mag_array[:, :, 1]

    # Downsample for quiver (every 2 cells)
    step = 2
    x_pos = np.arange(0, grid_spec.nx, step) * grid_spec.dx * 1e9 + grid_spec.dx * 1e9 / 2
    y_pos = np.arange(0, grid_spec.ny, step) * grid_spec.dy * 1e9 + grid_spec.dy * 1e9 / 2
    X, Y = np.meshgrid(x_pos, y_pos)

    U = mx[::step, ::step].T
    V = my[::step, ::step].T

    ax.quiver(X, Y, U, V, color='black', alpha=0.6, scale=15, width=0.003)

    plt.colorbar(im, ax=ax, label='mz')

plt.suptitle(f'FDM Simulation: mz Component Evolution (20x20x1 grid, dt={timeStep*1e12:.2f}ps)',
             fontsize=14, fontweight='bold')
plt.tight_layout()

# Save figure
output_file = 'fdm_2d_mz_snapshots.png'
plt.savefig(output_file, dpi=150, bbox_inches='tight')
print(f"\nPlot saved to: {output_file}")

# Create animation
print("\nGenerating animation...")
fig_anim, ax_anim = plt.subplots(figsize=(8, 8))

def update_frame(frame_idx):
    ax_anim.clear()
    mag_array = mags_to_array(snapshots[frame_idx], grid_spec.nx, grid_spec.ny)
    mz = mag_array[:, :, 2]
    mx = mag_array[:, :, 0]
    my = mag_array[:, :, 1]

    im = ax_anim.imshow(mz.T, cmap='RdBu_r', vmin=-1, vmax=1, origin='lower',
                        extent=[0, width*1e9, 0, length*1e9])

    # Add quiver
    step = 2
    x_pos = np.arange(0, grid_spec.nx, step) * grid_spec.dx * 1e9 + grid_spec.dx * 1e9 / 2
    y_pos = np.arange(0, grid_spec.ny, step) * grid_spec.dy * 1e9 + grid_spec.dy * 1e9 / 2
    X, Y = np.meshgrid(x_pos, y_pos)
    U = mx[::step, ::step].T
    V = my[::step, ::step].T
    ax_anim.quiver(X, Y, U, V, color='black', alpha=0.6, scale=15, width=0.003)

    ax_anim.set_title(f'FDM Magnetization Evolution: t = {snapshot_times[frame_idx]*1e9:.2f} ns')
    ax_anim.set_xlabel('x (nm)')
    ax_anim.set_ylabel('y (nm)')
    return [im]

anim = FuncAnimation(fig_anim, update_frame, frames=len(snapshots),
                     interval=100, blit=False, repeat=True)

# Save animation
anim_file = 'fdm_2d_mz_animation.gif'
anim.save(anim_file, writer='pillow', fps=10, dpi=100)
print(f"Animation saved to: {anim_file}")

plt.show()

print("\n" + "=" * 60)
print("Visualization complete!")
print("=" * 60)
