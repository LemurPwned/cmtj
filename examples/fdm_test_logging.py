"""
Test script for FDM magnetization logging feature.
Demonstrates how to log magnetization snapshots during simulation.
"""
import cmtj
import numpy as np
import matplotlib.pyplot as plt

# Create a simple FDM layer
layer_template = cmtj.Layer(
    id="free",
    mag=cmtj.CVector(0, 0, 1),  # Initial magnetization (normalized)
    anis=cmtj.CVector(0, 0, 1),  # Anisotropy axis
    Ms=1.4e6,  # Saturation magnetization (A/m)
    thickness=1e-9,  # Layer thickness (m)
    cellSurface=1e-18,  # Fake surface area
    demagTensor=[cmtj.CVector(0,0,0), cmtj.CVector(0,0,0), cmtj.CVector(0,0,0)],
    damping=0.01  # Gilbert damping
)
# Set anisotropy
K1 = 1.3e6  # Anisotropy constant (J/m^3)
layer_template.setAnisotropyDriver(cmtj.ScalarDriver.getConstantDriver(K1))

# Create FDM grid: 10nm x 10nm x 1nm with 1nm cell size
grid_spec = cmtj.fdm.FDMGridSpec(
    width=10e-9,
    length=10e-9,
    thickness=1e-9,
    cellSizeXY=1e-9,
    nz=1
)

print(f"Grid: {grid_spec.nx}x{grid_spec.ny}x{grid_spec.nz} = {grid_spec.nx * grid_spec.ny * grid_spec.nz} cells")

# Create FDM layer and junction
fdm_layer = cmtj.fdm.FDMLayer(layer_template, grid_spec)
fdm_junction = cmtj.fdm.FDMJunction([fdm_layer])

# Set initial magnetization pattern: vortex-like
nx, ny = grid_spec.nx, grid_spec.ny
init_mag = []
center_x, center_y = nx / 2.0, ny / 2.0
for iz in range(grid_spec.nz):
    for iy in range(ny):
        for ix in range(nx):
            dx = ix - center_x
            dy = iy - center_y
            r = np.sqrt(dx**2 + dy**2)
            if r < 0.1:
                # Center core pointing down
                m = cmtj.CVector(0, 0, -1)
            else:
                # Tangential circulation
                theta = np.arctan2(dy, dx)
                mx = -np.sin(theta)
                my = np.cos(theta)
                mz = 0.2 * np.exp(-r / 2)  # Slight z-component
                norm = np.sqrt(mx**2 + my**2 + mz**2)
                m = cmtj.CVector(mx/norm, my/norm, mz/norm)
            init_mag.append(m)

fdm_junction.setLayerMagnetisationGrid("free", init_mag)

# Set exchange interaction
A = 1e-11  # Exchange stiffness (J/m)
fdm_junction.setLayerExchangeStiffness("free", A)

# Apply external field pulse to excite dynamics
field_mag = 0.05 / (4 * np.pi * 1e-7 * layer_template.Ms)  # ~50 mT equivalent
field_driver = cmtj.AxialDriver(
    cmtj.NullDriver(), cmtj.NullDriver(),
    cmtj.ScalarDriver.getConstantDriver(field_mag)
)
fdm_junction.setLayerExternalFieldDriver("free", field_driver)

# Storage for logged data
logged_times = []
logged_mags = []

def log_callback(time, iteration, mags):
    """Callback function to log magnetization at each write frequency"""
    logged_times.append(time)
    # Convert C++ vector to numpy array for easier analysis
    layer_mags = np.array([[m.x, m.y, m.z] for m in mags[0]])
    logged_mags.append(layer_mags)
    print(f"  Logged t={time*1e9:.2f}ns (iteration {iteration})")

# Set the logging callback
fdm_junction.setLogCallback(log_callback)

# Run simulation
print("\nRunning simulation...")
totalTime = 1e-9  # 1 ns
timeStep = 1e-13  # 100 fs
writeFrequency = 1e-10  # 100 ps (10 snapshots)

fdm_junction.runSimulation(
    totalTime=totalTime,
    timeStep=timeStep,
    writeFrequency=writeFrequency,
    verbose=True,
    solverMode=cmtj.SolverMode.RK4
)

print(f"\nLogged {len(logged_times)} snapshots")
print(f"Times: {[f'{t*1e9:.1f}' for t in logged_times[:5]]} ... ns")

# Analyze magnetization dynamics
print("\nMagnetization statistics at each snapshot:")
for i, (t, mags) in enumerate(zip(logged_times, logged_mags)):
    avg_mz = np.mean(mags[:, 2])
    std_mz = np.std(mags[:, 2])
    print(f"  t={t*1e9:4.1f}ns: <mz>={avg_mz:+.3f} ± {std_mz:.3f}")

# Plot magnetization evolution
fig, axes = plt.subplots(2, 3, figsize=(12, 8))
fig.suptitle(f"FDM Magnetization Logging Demo (Exchange A={A:.1e} J/m)")

# Show snapshots at 6 different times
snapshot_indices = [0, 2, 4, 6, 8, len(logged_times)-1]
for idx, snap_idx in enumerate(snapshot_indices):
    ax = axes[idx // 3, idx % 3]
    if snap_idx < len(logged_mags):
        mags = logged_mags[snap_idx].reshape(ny, nx, 3)
        t = logged_times[snap_idx]

        # Plot mz component as color
        im = ax.imshow(mags[:, :, 2], cmap='RdBu_r', vmin=-1, vmax=1, origin='lower')
        ax.set_title(f't = {t*1e9:.2f} ns')
        ax.set_xlabel('x (nm)')
        ax.set_ylabel('y (nm)')

        # Overlay in-plane arrows
        skip = max(1, nx // 8)
        X, Y = np.meshgrid(range(nx), range(ny))
        U = mags[:, :, 0]
        V = mags[:, :, 1]
        ax.quiver(X[::skip, ::skip], Y[::skip, ::skip],
                 U[::skip, ::skip], V[::skip, ::skip],
                 color='black', scale=15, width=0.003, alpha=0.6)

        plt.colorbar(im, ax=ax, label='mz')

plt.tight_layout()
plt.savefig('/Users/jm/repos/cmtj/examples/fdm_logging_demo.png', dpi=150)
print(f"\nSaved plot to: examples/fdm_logging_demo.png")

# Clear callback (good practice)
fdm_junction.clearLogCallback()
print("\nDone!")
