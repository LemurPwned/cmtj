"""
FDM 2D Visualization - Fast Test Version
==========================================

Quick version for testing with coarser time step.
"""

import numpy as np
import matplotlib.pyplot as plt
import cmtj
from cmtj.utils import get_full_demag_tensor, convert_tensor_to_cpp
import time

# Physical parameters
Ms = 1e6
Ku = 1e4
damping = 0.05
thickness = 2e-9
A_exchange = 1e-11  # Exchange stiffness [J/m], typical for permalloy

# Grid: 20x20x1
width = 100e-9
length = 100e-9
cellSizeXY = 5e-9
nz = 1

# Fast simulation parameters
totalTime = 2e-9  # 2 ns
timeStep = 1e-13  # 0.1 ps (10x larger than requested, but faster)
n_segments = 8

print("="*60)
print("FDM 2D Visualization - Fast Test")
print("="*60)
print(f"Grid: 20x20x1, dt={timeStep*1e12}ps")
print(f"Total time: {totalTime*1e9}ns, segments: {n_segments}")
print(f"Total steps: {int(totalTime/timeStep):,}")
print("="*60)

# Create layer
layer_id = "free"
demagTensor = [cmtj.CVector(0,0,0), cmtj.CVector(0,0,0), cmtj.CVector(0,0,1)]

layer = cmtj.Layer(
    id=layer_id,
    mag=cmtj.CVector(0,0,1),
    anis=cmtj.CVector(0,0,1),
    Ms=Ms,
    thickness=thickness,
    cellSurface=cellSizeXY**2,
    demagTensor=demagTensor,
    damping=damping
)
layer.setAnisotropyDriver(cmtj.constantDriver(Ku))

# FDM grid
grid = cmtj.fdm.FDMGridSpec(width, length, thickness, cellSizeXY, nz)
fdm_layer = cmtj.fdm.FDMLayer(layer, grid)

# Set exchange stiffness
fdm_layer.setExchangeStiffness(A_exchange)
print(f"Exchange stiffness set: A = {A_exchange:.2e} J/m")

# Demag tensor
print("\nComputing demag tensor...")
t0 = time.time()
tensor = get_full_demag_tensor((grid.nx, grid.ny, grid.nz),
                                (grid.dx, grid.dy, grid.dz))
fdm_layer.setDemagTensor(convert_tensor_to_cpp(tensor))
print(f"Done in {time.time()-t0:.1f}s")

# Initial state: vortex
print("Setting vortex state...")
mags = []
cx, cy = grid.nx/2, grid.ny/2
for iz in range(grid.nz):
    for iy in range(grid.ny):
        for ix in range(grid.nx):
            dx, dy = ix-cx, iy-cy
            r = np.sqrt(dx**2 + dy**2)
            if r < 4:
                angle = np.arctan2(dy, dx)
                mx, my, mz = 0.5*np.cos(angle), 0.5*np.sin(angle), 0.7
            else:
                mx, my, mz = 0, 0, 1
            norm = np.sqrt(mx**2 + my**2 + mz**2)
            mags.append(cmtj.CVector(mx/norm, my/norm, mz/norm))
fdm_layer.setMagnetisationGrid(mags)

# Junction
junction = cmtj.fdm.FDMJunction([fdm_layer])

# External field
H = 1e5
junction.setLayerExternalFieldDriver(layer_id, cmtj.AxialDriver(
    cmtj.constantDriver(H), cmtj.NullDriver(), cmtj.constantDriver(H*0.1)
))

# Run simulation
print("\nRunning simulation...")
snapshots = [junction.getLayerMagnetisationGrid(layer_id)]
times = [0]

for i in range(n_segments):
    t0 = time.time()
    junction.runSimulation(totalTime/n_segments, timeStep,
                          totalTime/n_segments, False, cmtj.RK4)
    snapshots.append(junction.getLayerMagnetisationGrid(layer_id))
    times.append((i+1)*totalTime/n_segments)
    print(f"  Segment {i+1}/{n_segments} ({time.time()-t0:.1f}s, t={times[-1]*1e9:.2f}ns)")

print("Simulation done!")

# Plot
def to_array(mags, nx, ny):
    arr = np.zeros((nx, ny, 3))
    idx = 0
    for iy in range(ny):
        for ix in range(nx):
            arr[ix, iy] = [mags[idx].x, mags[idx].y, mags[idx].z]
            idx += 1
    return arr

print("\nGenerating plots...")
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
axes = axes.flatten()

plot_idx = [0, 2, 4, 5, 6, 8]
for i, idx in enumerate(plot_idx):
    if idx >= len(snapshots):
        idx = len(snapshots) - 1

    arr = to_array(snapshots[idx], grid.nx, grid.ny)
    mz = arr[:,:,2]

    ax = axes[i]
    im = ax.imshow(mz.T, cmap='RdBu_r', vmin=-1, vmax=1, origin='lower',
                   extent=[0, width*1e9, 0, length*1e9])

    # Quiver
    mx, my = arr[:,:,0], arr[:,:,1]
    step = 2
    xp = np.arange(0, grid.nx, step) * grid.dx * 1e9
    yp = np.arange(0, grid.ny, step) * grid.dy * 1e9
    X, Y = np.meshgrid(xp, yp)
    ax.quiver(X, Y, mx[::step,::step].T, my[::step,::step].T,
              color='white', alpha=0.8, scale=10)

    ax.set_title(f't={times[idx]*1e9:.2f}ns')
    ax.set_xlabel('x (nm)')
    ax.set_ylabel('y (nm)')
    plt.colorbar(im, ax=ax, label='$m_z$')

plt.suptitle(f'FDM: 20x20x1 grid, dt={timeStep*1e12}ps', fontsize=14)
plt.tight_layout()

plt.savefig('fdm_test.png', dpi=150)
print("Saved: fdm_test.png")
plt.close()

print("\nDone!")
