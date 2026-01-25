"""
FDM simulation: 30x30x1 grid for 20ns with magnetization logging.
Demonstrates vortex dynamics with exchange interaction.
"""
import cmtj
import numpy as np
import matplotlib.pyplot as plt
import time as pytime

print("=" * 70)
print("FDM Simulation: 30x30x1 grid, 20ns duration")
print("=" * 70)

# Create layer template
layer_template = cmtj.Layer(
    id="free",
    mag=cmtj.CVector(0, 0, 1),
    anis=cmtj.CVector(0, 0, 1),
    Ms=1.4e6,  # A/m
    thickness=1e-9,  # 1 nm
    cellSurface=1e-18,
    demagTensor=[cmtj.CVector(0,0,0), cmtj.CVector(0,0,0), cmtj.CVector(0,0,0)],
    damping=0.01
)
K1 = 1.3e6  # J/m^3
layer_template.setAnisotropyDriver(cmtj.ScalarDriver.getConstantDriver(K1))

# Create 30x30x1 grid with 1nm cells (30nm x 30nm x 1nm total)
grid_spec = cmtj.fdm.FDMGridSpec(
    width=30e-9,
    length=30e-9,
    thickness=1e-9,
    cellSizeXY=1e-9,
    nz=1
)

print(f"\nGrid configuration:")
print(f"  Dimensions: {grid_spec.nx}x{grid_spec.ny}x{grid_spec.nz}")
print(f"  Total cells: {grid_spec.nx * grid_spec.ny * grid_spec.nz}")
print(f"  Cell size: {grid_spec.dx*1e9:.2f}nm x {grid_spec.dy*1e9:.2f}nm x {grid_spec.dz*1e9:.2f}nm")
print(f"  Physical size: {grid_spec.width*1e9:.1f}nm x {grid_spec.length*1e9:.1f}nm x {grid_spec.thickness*1e9:.1f}nm")

# Create FDM layer and junction
fdm_layer = cmtj.fdm.FDMLayer(layer_template, grid_spec)
fdm_junction = cmtj.fdm.FDMJunction([fdm_layer])

# Initialize with vortex pattern
print("\nInitializing vortex magnetization pattern...")
nx, ny = grid_spec.nx, grid_spec.ny
init_mag = []
center_x, center_y = nx / 2.0, ny / 2.0

for iz in range(grid_spec.nz):
    for iy in range(ny):
        for ix in range(nx):
            dx = ix - center_x
            dy = iy - center_y
            r = np.sqrt(dx**2 + dy**2)
            
            if r < 1.5:
                # Core pointing down
                m = cmtj.CVector(0, 0, -1)
            else:
                # Tangential circulation with out-of-plane component
                theta = np.arctan2(dy, dx)
                mx = -np.sin(theta)
                my = np.cos(theta)
                mz = 0.3 * np.exp(-r / 5)
                norm = np.sqrt(mx**2 + my**2 + mz**2)
                m = cmtj.CVector(mx/norm, my/norm, mz/norm)
            init_mag.append(m)

fdm_junction.setLayerMagnetisationGrid("free", init_mag)

# Set exchange interaction
A = 1e-11  # J/m
fdm_junction.setLayerExchangeStiffness("free", A)
print(f"  Exchange stiffness: A = {A:.1e} J/m")

# Apply small external field to excite dynamics
field_mag = 0.02 / (4 * np.pi * 1e-7 * layer_template.Ms)
field_driver = cmtj.AxialDriver(
    cmtj.NullDriver(), cmtj.NullDriver(),
    cmtj.ScalarDriver.getConstantDriver(field_mag)
)
fdm_junction.setLayerExternalFieldDriver("free", field_driver)
print(f"  External field: ~20 mT in z-direction")

# Storage for logged data
logged_times = []
logged_mags = []
log_start_time = pytime.time()

def log_callback(time, iteration, mags):
    """Callback to log magnetization snapshots"""
    logged_times.append(time)
    layer_mags = np.array([[m.x, m.y, m.z] for m in mags[0]])
    logged_mags.append(layer_mags)
    elapsed = pytime.time() - log_start_time
    print(f"  t={time*1e9:5.1f}ns (iter {iteration:6d}) | Elapsed: {elapsed:.1f}s")

# Set logging callback
fdm_junction.setLogCallback(log_callback)

# Simulation parameters
totalTime = 20e-9  # 20 ns
timeStep = 1e-13   # 100 fs
writeFrequency = 1e-9  # 1 ns (20 snapshots)

print(f"\nSimulation parameters:")
print(f"  Total time: {totalTime*1e9:.0f} ns")
print(f"  Time step: {timeStep*1e12:.0f} fs")
print(f"  Write frequency: {writeFrequency*1e9:.0f} ns")
print(f"  Total iterations: {int(totalTime/timeStep):,}")
print(f"  Expected snapshots: {int(totalTime/writeFrequency) + 1}")

print("\nStarting simulation...")
print("-" * 70)
sim_start = pytime.time()

fdm_junction.runSimulation(
    totalTime=totalTime,
    timeStep=timeStep,
    writeFrequency=writeFrequency,
    verbose=True,
    solverMode=cmtj.SolverMode.RK4
)

sim_end = pytime.time()
print("-" * 70)
print(f"Simulation completed in {sim_end - sim_start:.1f}s")

# Analysis
print(f"\nLogged {len(logged_times)} snapshots")
print("\nMagnetization statistics:")
for i in range(0, len(logged_times), max(1, len(logged_times)//10)):
    t = logged_times[i]
    mags = logged_mags[i]
    avg_mx = np.mean(mags[:, 0])
    avg_my = np.mean(mags[:, 1])
    avg_mz = np.mean(mags[:, 2])
    std_mz = np.std(mags[:, 2])
    print(f"  t={t*1e9:5.1f}ns: <m> = ({avg_mx:+.3f}, {avg_my:+.3f}, {avg_mz:+.3f}) | σ_z={std_mz:.3f}")

# Create visualization
print("\nGenerating visualization...")
fig = plt.figure(figsize=(16, 10))
gs = fig.add_gridspec(3, 4, hspace=0.3, wspace=0.3)

# Plot snapshots at different times
snapshot_indices = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20]
snapshot_indices = [i for i in snapshot_indices if i < len(logged_mags)]

for idx, snap_idx in enumerate(snapshot_indices[:12]):
    ax = fig.add_subplot(gs[idx // 4, idx % 4])
    
    mags = logged_mags[snap_idx].reshape(ny, nx, 3)
    t = logged_times[snap_idx]
    
    # Plot mz component
    im = ax.imshow(mags[:, :, 2], cmap='RdBu_r', vmin=-1, vmax=1, 
                   origin='lower', extent=[0, nx, 0, ny])
    
    # Overlay in-plane vectors
    skip = max(1, nx // 10)
    X, Y = np.meshgrid(range(nx), range(ny))
    U = mags[:, :, 0]
    V = mags[:, :, 1]
    ax.quiver(X[::skip, ::skip], Y[::skip, ::skip],
             U[::skip, ::skip], V[::skip, ::skip],
             color='black', scale=20, width=0.003, alpha=0.7)
    
    ax.set_title(f't = {t*1e9:.1f} ns', fontsize=10)
    ax.set_xlabel('x (nm)')
    ax.set_ylabel('y (nm)')
    ax.set_aspect('equal')

# Add colorbar
fig.colorbar(im, ax=fig.get_axes(), label='mz', shrink=0.6, pad=0.02)

fig.suptitle(f'FDM Vortex Dynamics: {nx}x{ny}x{grid_spec.nz} grid, A={A:.1e} J/m, {totalTime*1e9:.0f}ns', 
             fontsize=14, fontweight='bold')

output_path = '/Users/jm/repos/cmtj/examples/fdm_30x30_20ns.png'
plt.savefig(output_path, dpi=150, bbox_inches='tight')
print(f"Saved visualization: {output_path}")

# Plot time evolution of average magnetization
fig2, axes = plt.subplots(1, 3, figsize=(15, 4))

times_ns = np.array(logged_times) * 1e9
mx_avg = [np.mean(m[:, 0]) for m in logged_mags]
my_avg = [np.mean(m[:, 1]) for m in logged_mags]
mz_avg = [np.mean(m[:, 2]) for m in logged_mags]

axes[0].plot(times_ns, mx_avg, 'b-', linewidth=2)
axes[0].set_xlabel('Time (ns)')
axes[0].set_ylabel('<mx>')
axes[0].grid(True, alpha=0.3)
axes[0].set_title('Average mx')

axes[1].plot(times_ns, my_avg, 'g-', linewidth=2)
axes[1].set_xlabel('Time (ns)')
axes[1].set_ylabel('<my>')
axes[1].grid(True, alpha=0.3)
axes[1].set_title('Average my')

axes[2].plot(times_ns, mz_avg, 'r-', linewidth=2)
axes[2].set_xlabel('Time (ns)')
axes[2].set_ylabel('<mz>')
axes[2].grid(True, alpha=0.3)
axes[2].set_title('Average mz')

plt.tight_layout()
output_path2 = '/Users/jm/repos/cmtj/examples/fdm_30x30_20ns_evolution.png'
plt.savefig(output_path2, dpi=150, bbox_inches='tight')
print(f"Saved evolution plot: {output_path2}")

# Cleanup
fdm_junction.clearLogCallback()
print("\n" + "=" * 70)
print("Simulation complete!")
print("=" * 70)
