# Solvers, temperature, and reproducibility

## SolverMode

`cmtj.SolverMode`: `RK4` (default), `DormandPrice` (adaptive step — note the API spelling, no "n"), `Heun`, `EulerHeun`.

```python
junction.runSimulation(5e-10, 1e-12, 1e-12, solverMode=SolverMode.EulerHeun)
```

- **RK4** — default, deterministic dynamics, fixed `dt`.
- **DormandPrice** — adaptive step, use instead of shrinking `dt` by hand for stiff/large-IEC systems.
- **Heun / EulerHeun — required for any stochastic run** (thermal noise, 1/f noise). RK4 cannot integrate the Stratonovich SDE correctly; if a temperature driver is set but the solver is left at RK4, results are wrong, not just noisier. `Heun` is the generally-preferred stochastic method (better convergence order than `EulerHeun`); see `docs/physics/macromagnetic_models.md` for why.

## Thermal / stochastic simulations

```python
junction.setLayerTemperatureDriver("free", constantDriver(300.0))  # K
junction.runSimulation(totalTime, dt, dt, solverMode=SolverMode.Heun)
```

- Needs a `ScalarDriver` on temperature (`constantDriver`, or time-varying for annealing protocols).
- Needs `Heun`/`EulerHeun` solver (see above).
- Stochastic noise magnitude scales with `1/sqrt(V)` (cell volume) — very small `cellSurface × thickness` volumes exaggerate thermal fluctuations; check the volume is physically sensible if switching looks too noisy or a p-bit never settles.

## Reproducibility — seeding

```python
junction.setLayerSeed("free", 1224)   # fixes the RNG for that layer
junction.setLayerSeed("free", None)   # re-randomizes it
```

Use this whenever a stochastic result needs to be reproducible (regression tests, curated examples, debugging a specific switching event). Without a seed, re-running the same script gives a different trajectory every time — expected for thermal noise, but a footgun if you're trying to compare "before/after" a code change and forgot to pin it. See `curated-examples/seeded-noise-demo.py` for the pattern (same seed twice vs. a different seed, plotted side by side).
