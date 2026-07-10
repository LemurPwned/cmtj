---
name: cmtj-simulation
description: Write correct CMTJ macromagnetic simulations (LLGS dynamics, Smit-Beljers statics, PIMM/VSD/FMR procedures, stochastic/thermal runs, parallel sweeps, stacks) with right units, layer setup, drivers, and solver choice. Use when writing/reviewing/debugging a cmtj simulation script, adding a Layer/Junction/Stack, setting anisotropy/IEC/field/SOT/STT/temperature drivers, choosing a SolverMode, parallelizing a sweep with cmtj.utils.parallel, or when values look off (switching won't trigger, frequencies wrong, unstable integration, stochastic results not reproducible).
---

# CMTJ simulation writing

CMTJ = macromagnetic sim for MTJs/spin Hall bars. C++ core (`core/`) + Python bindings (`_cmtj`), plus pure-Python static models (`cmtj/models/general_sb.py`).

**Read `AGENTS.md` (repo root) first** — full unit table, layer-naming rules, contribution/driver classification, PR conventions. This skill is the condensed, task-time checklist; AGENTS.md is the source of truth. If they ever disagree, AGENTS.md wins.

## Before writing anything

1. Find the closest analog in `curated-examples/` (hysteresis, switching, FMR, PIMM/VSD, IEC, p-bit, STO, AHE...) and copy its structure/parameter scale, don't invent from scratch.
2. Check `docs/physics/contributions.md` for the contribution you need (external field, Oersted, IEC, anisotropy, demag/dipole, thermal) — it has the exact setter name and unit.
3. Decide dynamic (LLGS, `core`/`Junction`) vs static (Smit-Beljers, `cmtj/models/general_sb.py`). Default to dynamic unless statics are explicitly requested.

## Units — check every value against this before running

| Quantity                  | Unit (core/Junction API)       | Unit (Smit-Beljers `general_sb.py`)            |
| ------------------------- | ------------------------------ | ---------------------------------------------- |
| `Ms`                      | **T** (0.5–1.6 typical)        | **A/m** — only exception, and spherical coords |
| `Ku`, anisotropy          | J/m³ (1e2–1e6)                 | J/m³                                           |
| `H_ext`, effective fields | A/m always                     | A/m                                            |
| `J` (IEC)                 | J/m² (±0.001–±3.0)             | J/m²                                           |
| `D` (iDMI)                | J/m² (0–3.0)                   | J/m²                                           |
| `thickness`               | m (0.8–2.0 **nm**, so `~1e-9`) | m                                              |

Common mistake: passing `Ms` in A/m to a `Layer` (should be T) — values ~1e6 instead of ~1. If a simulation looks physically wrong (no switching, wrong resonance frequency), check this first.

## Layer construction

- No SOT/STT → plain `Layer(...)`.
- STT/SOT needed → same `Layer`, but set `p` (polarisation) and STT/SOT drivers correctly.
- **Never create a second `Layer` for a fixed/reference layer.** Model it as a fixed polarisation direction on the free layer:
  ```python
  layer.setReferenceLayer(CVector(0, 1.0, 0.0))
  ```
- Naming: MTJ layers are `free` (tracked) and `reference`/`bottom`/`top` (fixed). Spin Hall bars: `topFM`/`bottomFM`, or `layerFM1`, `layerFM2`, ... for >2. Follow this even for internal variable names.

## Drivers — 3 categories (see `docs/physics/contributions.md`)

- `ScalarDriver`-driven (magnitude only, axis fixed elsewhere): anisotropy `H_K`, thermal, IEC.
- `AxialDriver`-driven (x/y/z magnitude independently): external field `H_ext`, Oersted `H_Oe`. Built from 3 `ScalarDriver`s.
- Static, not driver-based: demag/dipole tensors, passed once at `Layer` construction.

```python
from cmtj import AxialDriver, ScalarDriver, constantDriver

layer.setAnisotropyDriver(constantDriver(Ku))          # ScalarDriver contribution
junction.setLayerExternalFieldDriver(                    # AxialDriver contribution
    "free",
    AxialDriver(
        ScalarDriver.getConstantDriver(Hx),
        ScalarDriver.getConstantDriver(Hy),
        ScalarDriver.getConstantDriver(0),
    ),
)
```

Driver factory functions available beyond constant: `sineDriver`, `posSineDriver`, `pulseDriver`, `stepDriver`, `trapezoidDriver`, `gaussianImpulseDriver`, `gaussianStepDriver` — use these instead of hand-rolling time-dependent loops.

## Integration

- Fixed step `dt = 1e-12` s by default. Drop to `1e-13` or lower only for large IEC (`J` near or above 1e-4 J/m²) or if convergence is in question — confirm by halving `dt` and checking the result doesn't change.
- `SolverMode.RK4` is the default deterministic solver. `SolverMode.DormandPrice` (adaptive step — that's the real API spelling, no "n") only when explicitly warranted instead of hand-shrinking `dt`.
- **Any simulation with a temperature driver (`setLayerTemperatureDriver`) or noise must use `SolverMode.Heun` or `SolverMode.EulerHeun`.** RK4 does not correctly integrate the stochastic (Stratonovich) term — this produces silently wrong results, not just noisier ones. Details + reproducible seeding via `junction.setLayerSeed(...)`: `references/solvers-and-stochastic.md`.
- Sensible total sim time: 1–500 ns. Longer rarely makes physical sense here.
- Use `junction.clearLog()` / `stack.clearLogs()` before repeated runs (parameter scans) — keeps memory bounded and speeds up scans a lot.
- When scanning a parameter (field sweep, frequency sweep) **serially**, seed each step from the previous step's final state, not from scratch — faster convergence. Perturb slightly (don't reuse the exact previous vector) to avoid getting stuck in a local minimum, unless computing spectra only (where it doesn't matter). Running the sweep in parallel instead? See the relaxation caveat in `references/parallelism.md` — this trick does not carry over as-is.

## Don't hand-roll what the library already does

`cmtj.utils` has ready-made FFT/PIMM/VSD procedures, resistance formulas (AMR/SMR/AHE/GMR), and
field-vector helpers (`FieldScan`). A bespoke FFT loop or resistance formula is the most common
source of avoidable bugs (wrong windowing, wrong sign convention). Full list with import paths:
`references/procedures-and-utils.md`.

Sweeping 2+ independent parameters (e.g. field × frequency)? `cmtj.utils.parallel.distribute`
parallelizes it, but naively porting a serial sweep breaks the "seed from previous step" trick
above — see `references/parallelism.md` before parallelizing.

Doing statics with `cmtj/models/general_sb.py` (Smit-Beljers)? `Ms` is A/m and coordinates are
spherical there — the one unit exception in the whole codebase. Full pattern (`LayerSB`, `Solver`,
reusing the previous solve as the next initial guess): `references/smit-beljers.md`.

## Plotting (when the task wants a figure)

- Style with `scienceplots` (`with plt.style.context(["science", "no-latex"])`), guarded by `contextlib.suppress(ImportError)` on the import — see any `curated-examples/*.py`.
- Plot `M(H)` or `R(H)` with the component labeled (`m_x`, `R_xy`, ...), not an unlabeled generic curve.
- Trajectories go on a sphere: mostly transparent (low alpha), must stay readable.

## After writing

- Touched `core/` or `python/cmtj.cpp`? Run `pip install -e .` — editable installs don't recompile C++.
- New reusable simulation method/pattern? Add a curated example following the existing file header style (physics description docstring, sensible names, a result figure) — see `AGENTS.md` §5.
- Run `pytest` (and `tests/test_curated_examples.py` specifically if you touched a curated example; `UPDATE_REFERENCE=1 pytest tests/test_curated_examples.py` to regenerate reference data intentionally).
- Optional: `python skills/cmtj-simulation/scripts/check_units.py <file.py>` — heuristic static check that flags the most common unit/parameter mistakes (Ms magnitude, thickness scale, dt vs IEC magnitude, a layer named reference/fixed/pinned without `setReferenceLayer`). Not a substitute for reading the AGENTS.md rules, just a fast pre-flight grep. `scripts/post_edit_hook.py` wires the same check into a Claude Code `PostToolUse` hook (see its docstring for the `.claude/settings.json` snippet) if you want it to run automatically after every edit to a file that imports `cmtj`.
