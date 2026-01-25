"""
FDM (Finite Difference Method) Example
=======================================

This example demonstrates how to use the FDM module to simulate
a magnetic thin film with spatial discretization.

We'll simulate a simple ferromagnetic layer with:
- 3D grid discretization (10x10x1 cells)
- Applied external field
- Demagnetization field computed from demag tensors
"""

import numpy as np
import cmtj
from cmtj.utils import get_full_demag_tensor, convert_tensor_to_cpp

# Physical parameters
Ms = 1e6  # Saturation magnetization [A/m]
Ku = 1e3  # Anisotropy constant [J/m^3]
damping = 0.01  # Gilbert damping
thickness = 2e-9  # Layer thickness [m]

# Grid parameters
width = 50e-9  # Total width in x [m]
length = 50e-9  # Total length in y [m]
cellSizeXY = 5e-9  # Cell size in xy plane [m]
nz = 1  # Number of cells in z direction

# Create a base Layer template with physical parameters
layer_id = "free_layer"
initial_mag = cmtj.CVector(0, 0, 1)  # Initial magnetization along +z
anisotropy_axis = cmtj.CVector(0, 0, 1)  # Anisotropy along +z
cellSurface = cellSizeXY * cellSizeXY  # Will be overridden by FDMLayer

# DemagTensor for Layer (will be overridden by FDM layer's per-cell tensors)
# For now, provide a simple dipolar approximation
demagTensor = [
    cmtj.CVector(0, 0, 0),  # Nxx, Nxy, Nxz
    cmtj.CVector(0, 0, 0),  # Nyx, Nyy, Nyz
    cmtj.CVector(0, 0, 1)   # Nzx, Nzy, Nzz (simplified thin film approx)
]

layer_template = cmtj.Layer(
    id=layer_id,
    mag=initial_mag,
    anis=anisotropy_axis,
    Ms=Ms,
    thickness=thickness,
    cellSurface=cellSurface,
    demagTensor=demagTensor,
    damping=damping
)

# Set anisotropy constant
layer_template.setAnisotropyDriver(cmtj.constantDriver(Ku))

# Create FDM grid specification
grid_spec = cmtj.fdm.FDMGridSpec(
    width=width,
    length=length,
    thickness=thickness,
    cellSizeXY=cellSizeXY,
    nz=nz
)

print(f"Grid dimensions: {grid_spec.nx} x {grid_spec.ny} x {grid_spec.nz}")
print(f"Cell size: dx={grid_spec.dx*1e9:.2f}nm, dy={grid_spec.dy*1e9:.2f}nm, dz={grid_spec.dz*1e9:.2f}nm")
print(f"Total cells: {grid_spec.nx * grid_spec.ny * grid_spec.nz}")

# Create FDM layer
fdm_layer = cmtj.fdm.FDMLayer(layer_template, grid_spec)

# Compute and set demagnetization tensor
print("\nComputing demagnetization tensor...")
n = (grid_spec.nx, grid_spec.ny, grid_spec.nz)
dx = (grid_spec.dx, grid_spec.dy, grid_spec.dz)
demag_tensor = get_full_demag_tensor(n, dx)
print(f"Demag tensor shape: {demag_tensor.shape}")

# Convert to C++ format and set on layer
cpp_demag_tensor = convert_tensor_to_cpp(demag_tensor)
fdm_layer.setDemagTensor(cpp_demag_tensor)
print("Demagnetization tensor set successfully")

# Set initial magnetization with small perturbation
# Start mostly along +z with small tilt towards +x
mags = []
for iz in range(grid_spec.nz):
    for iy in range(grid_spec.ny):
        for ix in range(grid_spec.nx):
            # Small tilt from z towards x (10 degrees)
            theta = np.radians(10)
            mx = np.sin(theta)
            my = 0.0
            mz = np.cos(theta)
            mags.append(cmtj.CVector(mx, my, mz))

fdm_layer.setMagnetisationGrid(mags)
print(f"Set initial magnetization for {len(mags)} cells")

# Create FDM junction with the layer
layers = [fdm_layer]
junction = cmtj.fdm.FDMJunction(layers)

# Set external field driver (apply field along +x direction to switch magnetization)
H_ext = 2e5  # External field [A/m] - strong enough to overcome anisotropy
field_driver = cmtj.AxialDriver(
    cmtj.constantDriver(H_ext),  # Constant field in x
    cmtj.NullDriver(),  # No field in y
    cmtj.NullDriver()  # No field in z
)
junction.setLayerExternalFieldDriver(layer_id, field_driver)

# Run simulation
print("\nRunning FDM simulation...")
totalTime = 5e-9  # 5 nanoseconds
timeStep = 1e-13  # 0.1 picoseconds
writeFrequency = 1e-11  # Write every 10 ps

junction.runSimulation(
    totalTime=totalTime,
    timeStep=timeStep,
    writeFrequency=writeFrequency,
    verbose=True,
    solverMode=cmtj.RK4
)

# Get final magnetization state
final_mags = junction.getLayerMagnetisationGrid(layer_id)
print(f"\nSimulation complete!")
print(f"Final magnetization (first cell): ({final_mags[0].x:.3f}, {final_mags[0].y:.3f}, {final_mags[0].z:.3f})")

# Compute average magnetization
avg_mx = np.mean([m.x for m in final_mags])
avg_my = np.mean([m.y for m in final_mags])
avg_mz = np.mean([m.z for m in final_mags])
print(f"Average magnetization: ({avg_mx:.3f}, {avg_my:.3f}, {avg_mz:.3f})")
