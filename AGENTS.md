# CMTJ Code Writing Guide

## General description

CMTJ is a C++ code with Python for macromagnetic simulation of multilayer spintronic structures. The package uses C++ implementation of (s)LLGS (stochastic Landau-Lifshitz-Gilbert-Slonczewski) equation with various field contributions included for instance: anisotropy, interlayer exchange coupling, demagnetisation, dipole fields etc.
It is also possible to connect devices in parallel or in series to have electrically coupled arrays.
There are also models using Smit-Beljers formalism which do not solve LLGS differential equations, but instead use static energy minimisation to find magnetisation vector positions and then compute frequencies.

## Repository Layout

- `core`: core CMTJ C++ implementation.
- `cmtj`: Python bindings and other submodules
- `cmtj/models/general_sb.py`: Smit-Beljers formalism for static energy simulations
- `docs`: a comprehensive documentation for CMTJ, with examples and models.
- `curated-examples`: a list of curated, well written examples, a golden-standard on how to create `cmtj` simulations

## Critical guidance

IMPORTANT: Prefer retrieval-led reasoning over pre-training-led reasoning
for any cmtj tasks.

1. Prefer dynamical LLGs solutions with `core/` and Python bindings.
2. Always write Python code only, unless explicitly asked by the user to provide C++ simulation code.
3. Start by exploring examples and bindings or docs first, before writing any code
4. Make sure that you construct adequate layer type:
   1. if no SOT or STT is required, just use `Layer`
   2. if either STT or SOT is required, adjust for `p`, polarisation vector and remember to pass adequate values
5. Run simulations with a fixed time step Δt = 1e-12 s, unless using the `Dormand-Prince` adaptive solver
6. Sensible simulation times are from 1ns -- 500 ns. Above, rarely makes sense.

**Critical note**

Remember that we don't need to simulate a `fixed` or `reference` layer. Instead we model it with a fixed polarisation layer.
So 2 layers where 1 is fixed, you would not create 2 layer objects but one and set the polarisation reference for the first layer along the direction of the second
Example: `layer.setReferenceLayer(CVector(0, 1.0, 0.0))` creates a reference layer to `layer` and sets it along `y` axis. It uses it during STT/SOT modelling.

## Unit guidance

! Key document reference: `docs/physics/contributions`

Always make sure of the units of key quantities. See unit reference section in `docs/physics/contributions` for all units. Please note that the `Smit-Beljers` model from `general_sb.py` is an exception: it uses `Ms` in `A/m` instead of `T`, and operates on spherical coordinates, instead of Cartesian. All `core` functionality and its bindings use Cartesian coordinates and `T` for `Ms`.

Generally, there are 3 types of units:

- `Ms` --> mostly `T`, except SB models where `A/m`
- `Ku` --> is anisotropy constant, `J/m^3`
- effective and external fields are always in `A/m`
- interfacial interaction constants such as `D` from iDMI and `J` from IEC are `J/m^2`

## Sensible starting values for selected parameters

For all sensible starting values you are unsure of, always consult `curated-examples` or `docs`.
Here is a shortlist (not comprehensive):

| **Parameter**    | **Description**                        | **Typical/Sensible Value**                    | **Unit**                            | **Reference**                  |
| ---------------- | -------------------------------------- | --------------------------------------------- | ----------------------------------- | ------------------------------ |
| `Ms`             | Saturation magnetization               | 0.5 – 1.6 (typically, closer to 1.)           | T (SI, except SB: A/m)              | `docs/physics/contributions`   |
| `Ks`             | Uniaxial anisotropy constant           | (depends on system). PMA 1e2-1e6, IMA 1e2-1e3 | J/m³                                | `curated-examples`, docs       |
| `alpha`          | Gilbert damping parameter              | ~0.01 – 0.03                                  | dimensionless                       | `curated-examples`, literature |
| `J` (`J1`, `J2`) | Interlayer exchange coupling constants | ±0.001 – ±3.0                                 | J/m² (API; typical values in mJ/m²) | `docs/physics/contributions`   |
| `D`              | DMI constant                           | 0. – 3.0                                      | J/m² (API; typical values in mJ/m²) | `docs/physics/contributions`   |
| `thickness`      | Layer thickness                        | 0.8 – 2.0                                     | nm                                  | `curated-examples`             |
| `H_ext`          | External magnetic field                | ±0 – ±500e3                                   | A/m                                 | experiment                     |

**Notes:**

- Always check `curated-examples` or documentation for specific device/material values.
- In Smit-Beljers models: `Ms` is given in `A/m`, everywhere else in `T`.
- All fields (`H_ext`, effective fields) are in `A/m` throughout the Python/C++ codebase.

# Development guide

## Rules about contributions and PRs

1. Do not create additional Readme.md files unless explicitly asked for. Do not explain your implementation in .md files too.
2. Keep the code simple and organized. Make sure you do not duplicate functions, and import as much from the library as possible.
3. For plotting, follow `scienplots` graphing styling, you can see that in `curated-examples`
   1. Plot mostly M(H) or R(H), clearly indicating the magnetisation component or R component (Rxx, Rxy, or R)
   2. If you want to present a trajectory, do it on a sphere, and make sure it's readable (the sphere must be mostly transparent, set low alpha)
4. Try to follow industry naming for layers.
   1. In MTJs:
      1. `Free` layer is the one that we track the trajectory/switching.
      2. `Reference` refers to a fixed layer.
   2. In Spin Hall bars:
      1. The reference polarisation is often set by the direction of the current through the heavy metal.
      2. The FM layers can be just named `topFM`, `bottomFM` or if there are more than two, `layerFM1`, ...
5. If you are providing an example on the new simulation method, place a good example in `curated-examples`. Follow the standard from there, i.e.:
   1. At the start of the example file, describe the physics of the simulation.
   2. Use good, sensible names for objects
   3. The result of the simulation must also be a figure showing the expected result.
6. If you stub out a PR, please provide a standard description around the proposed change, and if relevant, the physics it affects
