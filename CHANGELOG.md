# Changelog

## 1.11.0 

- Run of the mill optmisations + benchmarks. Even up to 40% speedup (best results with maximum optmisation flags, see README.md)

## 1.10.0 

- Added generalised constants which an be modified via the interface
- BREAKING: _dropped Python3.9-3.10_ support

## 1.9.1

- Switched to more modern builds with multilinux

## 1.9.0

- Added generalised linearisation options for frequency computing in linear models
- Revamped SB modelling for better processing and analytical roots + code cleanup

## 1.8.0

- Fixing bugs in electrical coupling
- You can now have useKCL = True Stacks which will implicitly preserve Kirchoff's current law, and will simulate perfect voltage sources

## 1.7.0

- Adaptive simulation is now available again with Dormand-Price method
- Adding curated simulation examples
- Memory cells are now supported (you can specify two polarisation vectors per Layer now)
- Import fixes and other minor bugfixes

## 1.6.2-1.6.6

- Fixed a bug in the `ScalarDriver` class which prevented the use of custom drivers with no arguments.

## 1.6.1

- Small fixes to the noise interfaces and resistance functions.
- Documentation updates: added a tutorial on custom dipole interactions + interface fixes and typing.

## 1.6.0

- Extended the `Stack` models allowing for non-symmetric coupling between devices.
  `Stack` current drivers can now be of any type and are adequately scaled.
- Custom definition of the `ScalarDriver` is now possible and documented.
- Fixed a bug in the `Stack` class which inverted the connection order of in-series connections.
- Exposed IDMI interaction to Layer and Junction classes.
- Added `getLayer` method to the `Junction` class and `getJunction` method to the `Stack` class that return a reference to the object.
- Fixed and expanded the `reservoir` module. Now, `GroupInteraction` can use any dipole interaction function, with 3 provided as default: `computeDipoleInteraction`, `computeDipoleInteractionNoumra` and `nullDipoleInteraction` (0 dipole tensor).

## 1.5.0-1.5.4

- Dipole interaction added to the `SB Model`
- Kasdin 1/f noise generator added to the `noise` module and to the solvers
- Reworking the solvers for better performance and stability
- Added a simple noise model to the `utils` class. It exists outside standard simulation procedures.
- Added LLGB bindings and code. The solver is still WIP and doesn't integrate with more advanced features yet.
- Added aliases for `ScalarDriver` -- for example, instead of calling `ScalarDriver.getConstantDriver`, you can now call `constantDriver` directly to create a constant driver.
- Improve stub detection across editors and IDEs

## 1.4.1

- Adding a basic optimisation script in the `optimization` module.
- Streamlit optimization updates.

## 1.4.0

- Adding new, dynamic symbolic model compatible with `Solver` class. It is now possible to use the `Solver` class with `LayerDynamic` to solve the LLG equation.
- Added tests for procedures and operators.
- Added missing operators for `CVector` class in Python.
- `CVector` is now subscriptable.
- Added new `CVector` and `AxialDriver` initialisations.
- `VSD` and `PIMM` procedures accept additional new parameters.
- Added some optimization utilities like `coordinate_descent`.
- Added a `streamlit` service for an example PIMM simulation.

## 1.3.2

- Added new `ScalarDrivers` -- Gaussian impulse and Gaussian step.

## 1.3.1

- `CVector` got extra functionality in Python bindings. Operators are now supported.
- Domain Wall dynamics is now also for 2 layer systems. Added edge potential.
- SB model generalised for N layers.

## 1.3.0

- Adding DW dynamics 1D model with dynamic drivers. (Numba optimised)
- Adding SB model for energy-based FMR computation. Gradient computed using Adam optimiser.
- Moving resistance functions from `utils` to `resistance`
- Introducing docs updates for tutorial notebook (dark/light toggle works now).
- Reservoir computing is now exposed in Python in the `reservoir` computing module.

## 1.2.0

- Oersted field computation helper class in [cmtj/models/oersted.py](cmtj/models/oersted.py). Basic functionality is there, but needs to be further tested and documented. Next release potentially will move the computation to C++ for speed.
- Added Heun (2nd order) solver and made it default for thermal computation. This is a more stable solver than the Euler solver, but is slower. The Euler solver is still available as an option.
- Stack class now supports arbitrary layer ids to be coupled.
- Extended the plotting capabilities of the Stack class. Now supports plotting of the magnetic field and the current density.
- Added alternative STT formulation which in some cases may be useful.
- Fixed some minor bugs in the thermal solver.
- Fixed some minor bugs in the Stack class.
- Updating tutorials on the docs page.
- Bunch of extra documentation and examples.
