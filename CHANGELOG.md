# Changelog

# 1.3.3

- Added tests for procedures.
- `VSD` and `PIMM` procedures accept additional new parameters.

# 1.3.2

- Added new `ScalarDrivers` -- Gaussian impulse and Gaussian step.

# 1.3.1

- `CVector` got extra functionality in Python bindings. Operators are now supported.
- Domain Wall dynamics is now also for 2 layer systems. Added edge potential.
- SB model generalised for N layers.

## 1.3.0

- Adding DW dynamics 1D model with dynamic drivers. (Numba optimised)
- Adding SB model for energy-based FMR computation. Gradient computed using Adam optimiser.
- Moving resistance functions from `utils` to `resistance`
- Introducting docs updates for tutorial notebook (dark/light toggle works now).
- Reservoir computing is now exposed in Python in the `reservoir` computing module.

## 1.2.0

- Oersted field computation helper class in [cmtj/models/oersted.py](cmtj/models/oersted.py). Basic functionality is there, but needs to be futher tested and documented. Next release potentially will move the computation to C++ for speed.
- Added Heun (2nd order) solver and made it default for thermal computation. This is a more stable solver than the Euler solver, but is slower. The Euler solver is still available as an option.
- Stack class now supports arbitrary layer ids to be coupled.
- Extended the plotting capabilities of the Stack class. Now supports plotting of the magnetic field and the current density.
- Added alternative STT formulation which in some cases may be useful.
- Fixed some minor bugs in the thermal solver.
- Fixed some minor bugs in the Stack class.
- Updating tutorials on the docs page.
- Bunch of extra documentation and examples.
