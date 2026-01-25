"""
FDM simulation: 30x30x1 grid for 20ns with STRONG excitation.
Demonstrates dynamic vortex motion with exchange interaction.
"""
import cmtj
import numpy as np
import matplotlib.pyplot as plt
import time as pytime

print("=" * 70)
print("FDM Vortex Dynamics: 30x30x1 grid, 20ns with STRONG excitation")
print("=" * 70)

# Create layer template with LOWER damping for more dynamics
layer_template = cmtj.Layer(
    id="free",
    mag=cmtj.CVector(0, 0, 1),
    anis=cmtj.CVector(0, 0, 1),
    Ms=1.4e6,  # A/m
    thickness=1e-9,  # 1 nm
    cellSurface=1e-18,
    demagTensor=[cmtj.CVector(0,0,0), cmtj.CVector(0,0,0), cmtj.CVector(0,0,0)],
    damping=0.005  # Lower damping = more dynamics
)
K1 = 5e5  # Reduced anisotropy J/m^3
layer_template.setAnisotropyDriver(cmtj.ScalarDriver.getConstantDriver(K1))

# Create 30x30x1 grid
grid_spec = cmtj.fdm.FDMGridSpec(
    width=30e-9,
    length=30e-9,
    thickness=1e-9,
    cellSizeXY=1e-9,
    nz=1
)

print(f"\nGrid: {grid_spec.nx}x{grid_spec.ny}x{grid_spec.nz} = {grid_spec.nx*grid_spec.ny*grid_spec.nz} cells")
print(f"Physical size: {grid_spec.width*1e9:.0f}nm x {grid_spec.length*1e9:.0f}nm x {grid_spec.thickness*1e9:.0f}nm")

# Create FDM layer and junction
fdm_layer = cmtj.fdm.FDMLayer(layer_template, grid_spec)
fdm_junction = cmtj.fdm.FDMJunction([fdm_layer])

# Initialize with OFFSET vortex to create dynamics
print("\nInitializing OFFSET vortex pattern...")
nx, ny = grid_spec.nx, grid_spec.ny
init_mag = []
center_x, center_y = nx / 2.0 + 5, ny / 2.0 - 3  # Offset from center

for iz in range(grid_spec.nz):
    for iy in range(ny):
        for ix in range(nx):
            dx = ix - center_x
            dy = iy - center_y
            r = np.sqrt(dx**2 + dy**2)
            
            if r < 2.0:
                # Core pointing down
                m = cmtj.CVector(0, 0, -1)
            else:
                # Tangential circulation
                theta = np.arctan2(dy, dx)
                mx = -np.sin(theta) * (1 - np.exp(-r/10))
                my = np.cos(theta) * (1 - np.exp(-r/10))
                mz = 0.5 * np.exp(-r / 4)
                norm = np.sqrt(mx**2 + my**2 + mz**2)
                m = cmtj.CVector(mx/norm, my/norm, mz/norm)
            init_mag.append(m)

fdm_junction.setLayerMagnetisationGrid("free", init_mag)

# Set exchange interaction
A = 1.3e-11  # Slightly higher J/m
fdm_junction.setLayerExchangeStiffness("free", A)
print(f"Exchange: A = {A:.1e} J/m")
print(f"Anisotropy: K1 = {K1:.1e} J/m^3")
print(f"Damping: α = {layer_template.damping}")

# Apply PULSED external field
field_mag = 0.15 / (4 * np.pi * 1e-7 * layer_template.Ms)  # ~150 mT
field_driver = cmtj.AxialDriver(
    cmtj.ScalarDriver.getPulseDriver(0, 2e-9, field_mag, field_mag),  # 2ns pulse in x
    cmtj.ScalarDriver.getPulseDriver(0, 2e-9, field_mag*0.5, field_mag*0.5),  # 2ns pulse in y
    cmtj.NullDriver()
)
fdm_junction.setLayerExternalFieldDriver("free", field_driver)
print(f"Field pulse: ~150 mT in xy-plane for 2ns")

# Logging
logged_times = []
logged_mags = []
vortex_centers = []  # Track vortex core position

def find_vortex_center(mags_2d):
    """Find vortex core as minimum mz location"""
    mz = mags_2d[:, :, 2]
    min_idx = np.argmin(mz)
    cy, cx = np.unravel_index(min_idx, mz.shape)
    return cx, cy

def log_callback(time, iteration, mags):
    logged_times.append(time)
    layer_mags = np.array([[m.x, m.y, m.z] for m in mags[0]])
    logged_mags.append(layer_mags)
    
    # Track vortex core
    mags_2d = layer_mags.reshape(ny, nx, 3)
    cx, cy = find_vortex_center(mags_2d)
    vortex_centers.append((cx, cy))
    
    if len(logged_times) % 5 == 0 or len(logged_times) <= 2:
        print(f"  t={time*1e9:5.1f}ns | Core at ({cx:.1f}, {cy:.1f})")

fdm_junction.setLogCallback(log_callback)

# Run simulation
totalTime = 20e-9
timeStep = 1e-13
writeFrequency = 0.5e-9  # 0.5ns intervals (40 snapshots)

print(f"\nSimulation: {totalTime*1e9:.0f}ns, dt={timeStep*1e12:.0f}fs, log every {writeFrequency*1e9:.1f}ns")
print(f"Total iterations: {int(totalTime/timeStep):,}")
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
print(f"Performance: {(sim_end - sim_start) / (totalTime*1e9):.2f} s per ns simulated")

# Analysis
print(f"\nLogged {len(logged_times)} snapshots")

# Vortex trajectory analysis
vortex_centers = np.array(vortex_centers)
displacements = np.sqrt(np.sum((vortex_centers - vortex_centers[0])**2, axis=1))
print(f"\nVortex core motion:")
print(f"  Initial position: ({vortex_centers[0,0]:.1f}, {vortex_centers[0,1]:.1f})")
print(f"  Final position: ({vortex_centers[-1,0]:.1f}, {vortex_centers[-1,1]:.1f})")
print(f"  Max displacement: {np.max(displacements):.1f} cells ({np.max(displacements)*1e-9*1e9:.1f} nm)")

# Create comprehensive visualization
print("\nGenerating plots...")
fig = plt.figure(figsize=(20, 12))
gs = fig.add_gridspec(4, 5, hspace=0.35, wspace=0.35)

# Plot snapshots
snapshot_indices = np.linspace(0, len(logged_mags)-1, 15, dtype=int)

for idx, snap_idx in enumerate(snapshot_indices):
    if idx >= 15:
        break
    ax = fig.add_subplot(gs[idx // 5, idx % 5])
    
    mags = logged_mags[snap_idx].reshape(ny, nx, 3)
    t = logged_times[snap_idx]
    cx, cy = vortex_centers[snap_idx]
    
    # Plot mz component
    im = ax.imshow(mags[:, :, 2], cmap='RdBu_r', vmin=-1, vmax=1, 
                   origin='lower', extent=[0, nx, 0, ny])
    
    # Mark vortex core
    ax.plot(cx, cy, 'ko', markersize=8, markerfacecolor='yellow', 
            markeredgewidth=2, label='Core')
    
    # In-plane vectors
    skip = 3
    X, Y = np.meshgrid(range(nx), range(ny))
    U = mags[:, :, 0]
    V = mags[:, :, 1]
    ax.quiver(X[::skip, ::skip], Y[::skip, ::skip],
             U[::skip, ::skip], V[::skip, ::skip],
             color='black', scale=25, width=0.002, alpha=0.6)
    
    ax.set_title(f't = {t*1e9:.1f} ns', fontsize=9, fontweight='bold')
    ax.set_xlabel('x (nm)', fontsize=8)
    ax.set_ylabel('y (nm)', fontsize=8)
    ax.tick_params(labelsize=7)

# Vortex trajectory plot
ax_traj = fig.add_subplot(gs[3, :2])
ax_traj.plot(vortex_centers[:, 0], vortex_centers[:, 1], 'b-', linewidth=2, alpha=0.6)
ax_traj.plot(vortex_centers[0, 0], vortex_centers[0, 1], 'go', markersize=12, 
             label='Start', zorder=5)
ax_traj.plot(vortex_centers[-1, 0], vortex_centers[-1, 1], 'ro', markersize=12, 
             label='End', zorder=5)
# Add time markers
time_markers = np.linspace(0, len(vortex_centers)-1, 5, dtype=int)
for i in time_markers[1:-1]:
    ax_traj.plot(vortex_centers[i, 0], vortex_centers[i, 1], 'ko', markersize=6, alpha=0.5)
ax_traj.set_xlabel('x (cells)', fontweight='bold')
ax_traj.set_ylabel('y (cells)', fontweight='bold')
ax_traj.set_title('Vortex Core Trajectory', fontweight='bold')
ax_traj.legend()
ax_traj.grid(True, alpha=0.3)
ax_traj.set_aspect('equal')

# Average magnetization evolution
ax_mag = fig.add_subplot(gs[3, 2:])
times_ns = np.array(logged_times) * 1e9
mx_avg = [np.mean(m[:, 0]) for m in logged_mags]
my_avg = [np.mean(m[:, 1]) for m in logged_mags]
mz_avg = [np.mean(m[:, 2]) for m in logged_mags]

ax_mag.plot(times_ns, mx_avg, 'b-', linewidth=2, label='<mx>', alpha=0.7)
ax_mag.plot(times_ns, my_avg, 'g-', linewidth=2, label='<my>', alpha=0.7)
ax_mag.plot(times_ns, mz_avg, 'r-', linewidth=2, label='<mz>', alpha=0.7)
ax_mag.axvspan(0, 2, alpha=0.2, color='yellow', label='Field pulse')
ax_mag.set_xlabel('Time (ns)', fontweight='bold')
ax_mag.set_ylabel('Average Magnetization', fontweight='bold')
ax_mag.set_title('Magnetization Evolution', fontweight='bold')
ax_mag.legend(loc='best')
ax_mag.grid(True, alpha=0.3)

fig.colorbar(im, ax=fig.get_axes()[:15], label='mz', shrink=0.8, pad=0.01)

fig.suptitle(f'FDM Vortex Dynamics: {nx}x{ny}x{grid_spec.nz} grid, A={A:.1e} J/m, α={layer_template.damping}, {totalTime*1e9:.0f}ns', 
             fontsize=16, fontweight='bold')

output_path = '/Users/jm/repos/cmtj/examples/fdm_30x30_20ns_dynamic.png'
plt.savefig(output_path, dpi=150, bbox_inches='tight')
print(f"\nSaved: {output_path}")

# Create animation-style frames if needed
print("\nCreating detailed time evolution plot...")
fig2, axes = plt.subplots(2, 2, figsize=(14, 12))

# Vortex core position vs time
ax = axes[0, 0]
ax.plot(times_ns, vortex_centers[:, 0], 'b-', linewidth=2, label='x position')
ax.plot(times_ns, vortex_centers[:, 1], 'r-', linewidth=2, label='y position')
ax.axvspan(0, 2, alpha=0.2, color='yellow', label='Field pulse')
ax.set_xlabel('Time (ns)', fontweight='bold')
ax.set_ylabel('Position (cells)', fontweight='bold')
ax.set_title('Vortex Core Position vs Time', fontweight='bold')
ax.legend()
ax.grid(True, alpha=0.3)

# Displacement from initial
ax = axes[0, 1]
ax.plot(times_ns, displacements, 'purple', linewidth=2)
ax.axvspan(0, 2, alpha=0.2, color='yellow', label='Field pulse')
ax.set_xlabel('Time (ns)', fontweight='bold')
ax.set_ylabel('Displacement (cells)', fontweight='bold')
ax.set_title('Core Displacement from Initial Position', fontweight='bold')
ax.grid(True, alpha=0.3)

# Magnetization components
ax = axes[1, 0]
ax.plot(times_ns, mx_avg, 'b-', linewidth=2, label='<mx>')
ax.plot(times_ns, my_avg, 'g-', linewidth=2, label='<my>')
ax.plot(times_ns, mz_avg, 'r-', linewidth=2, label='<mz>')
ax.axvspan(0, 2, alpha=0.2, color='yellow', label='Field pulse')
ax.set_xlabel('Time (ns)', fontweight='bold')
ax.set_ylabel('Magnetization', fontweight='bold')
ax.set_title('Average Magnetization Components', fontweight='bold')
ax.legend()
ax.grid(True, alpha=0.3)

# Final state
ax = axes[1, 1]
final_mags = logged_mags[-1].reshape(ny, nx, 3)
im = ax.imshow(final_mags[:, :, 2], cmap='RdBu_r', vmin=-1, vmax=1, origin='lower')
skip = 3
X, Y = np.meshgrid(range(nx), range(ny))
ax.quiver(X[::skip, ::skip], Y[::skip, ::skip],
         final_mags[::skip, ::skip, 0], final_mags[::skip, ::skip, 1],
         color='black', scale=25, width=0.003)
ax.plot(vortex_centers[-1, 0], vortex_centers[-1, 1], 'yo', markersize=10, 
        markeredgecolor='black', markeredgewidth=2)
ax.set_title(f'Final State (t={times_ns[-1]:.1f}ns)', fontweight='bold')
ax.set_xlabel('x (nm)')
ax.set_ylabel('y (nm)')
plt.colorbar(im, ax=ax, label='mz')

plt.tight_layout()
output_path2 = '/Users/jm/repos/cmtj/examples/fdm_30x30_20ns_analysis.png'
plt.savefig(output_path2, dpi=150, bbox_inches='tight')
print(f"Saved: {output_path2}")

# Cleanup
fdm_junction.clearLogCallback()

print("\n" + "=" * 70)
print("COMPLETE! Generated dynamic vortex evolution visualization.")
print("=" * 70)
