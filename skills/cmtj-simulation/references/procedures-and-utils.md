# Use the library helpers — don't hand-roll them

`cmtj.utils` has ready-made implementations for the measurement procedures and post-processing
that show up in almost every simulation. Writing a bespoke FFT loop or resistance formula when
one of these already exists is the most common source of avoidable bugs (wrong FFT windowing,
wrong AMR/SMR/AHE sign convention, etc.) — check here before writing one from scratch.

## Field vectors — `cmtj.utils.linear.FieldScan`

```python
from cmtj.utils import FieldScan

h = FieldScan.angle2vector(theta=89, phi=0.1, amplitude=100e3)   # -> CVector, angles in DEGREES
Hvecs = FieldScan.amplitude_scan(start, stop, steps, theta, phi)  # array of CVectors along one direction
thetas = FieldScan.theta_scan(start, stop, steps, amplitude, phi)
phis = FieldScan.phi_scan(start, stop, steps, amplitude, theta)
```

Note: `FieldScan`/Smit-Beljers angle helpers take **degrees**; the C++ core `AxialDriver` /
`CVector` API takes plain field components in A/m (no angle convention) — don't mix them up.

## PIMM / VSD measurement procedures — `cmtj.utils.procedures`

```python
from cmtj.utils.procedures import PIMM_procedure, VSD_procedure, ResistanceParameters

spectrum, freqs, other = PIMM_procedure(
    junction, Hvecs, int_step=1e-12,
    resistance_params=[ResistanceParameters(Rxx0=100, Rahe=1, w=..., l=...)],
)
```

- `PIMM_procedure` runs the pulsed-Oersted-field excitation + FFT for you (Pulse Induced Microwave
  Magnetometry), including the resistance log — don't reimplement the excitation-pulse +
  `compute_spectrum_strip` loop by hand.
- `VSD_procedure` is the equivalent for Voltage Spin Diode (sine-driven, mixed-down DC signal via
  `compute_sd`).
- `ResistanceParameters` is a plain dataclass (`Rxx0`, `Rxy0`, `Rahe`, `Rsmr`, `Ramr`, `w`, `l`) — fill
  in only the fields relevant to the device geometry, defaults are 0.

## Resistance — `cmtj.utils.resistance`

- `calculate_resistance_series` / `calculate_resistance_parallel` — GMR/TMR-style resistance from
  layer magnetizations for series- or parallel-connected stacks. Junction's own `R_<a>_<b>` log key
  (see `docs/physics/paralellism.md` example: `log['R_free_bottom']`) is produced by whichever
  `resistance_fn` you pass to a procedure.
- `calculate_magnetoresistance` / `compute_gmr` — plain `Rp`, `Rap` two-state MR, when you don't need
  the angular AMR/SMR/AHE decomposition.
- `Rxx_symbolic` / `Rxy_symbolic` / `calculate_linearised_resistance*` — symbolic/linearised forms for
  small-angle (harmonic Hall) analysis, used by `curated-examples/harmonic-hall.py` and `ahe_loops.py`.
- `compute_sd` — VSD DC mixing signal from a dynamic resistance trace and a reference current; used
  internally by `VSD_procedure`, but callable directly if you're not going through the procedure.

## Filters — `cmtj.utils.filters`

Post-processing on spectra (detrending, smoothing) before peak-fitting or plotting — see
`docs/tipsandtricks.md`: a plain `log`/`fft` is often noisy enough that peak-finding fails without
this step.
