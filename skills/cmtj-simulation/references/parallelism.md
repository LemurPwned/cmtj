# Parallel sweeps — `cmtj.utils.parallel.distribute`

Full walkthrough: `docs/physics/paralellism.md`.

## When to reach for this

A sweep with 2+ independent axes (e.g. frequency × field, or a parameter study across `J` values)
naively costs `O(N*K)`. `distribute` runs the outer axis across a worker pool so it collapses to
`O(K)` wall-clock.

```python
from cmtj.utils.parallel import distribute

def simulate_vsd(J, H, frequency):
    ...  # build junction, run, return a scalar/array result

Js = np.linspace(1e-4, 1e-3, 10)
Hrange = np.linspace(-15e3, 15e3, 100)
fscan = np.arange(1e9, 6.2e9, 0.2e9)
VSD = np.zeros((len(Js), len(fscan), len(Hrange)), dtype=np.float32)
for (k, i, j), out in distribute(simulate_vsd, [Js, Hrange, fscan]):
    VSD[k, j, i] = out
```

`distribute` yields `(indices, result)` pairs; `indices` match the order of the argument lists you
passed in (here: `J`, `H`, `frequency`), not the order of the output array — index carefully.

## The caveat that bites people

In a **serial** field/frequency scan, you get a "free" speedup by seeding each step from the
previous step's relaxed magnetization state (see main `SKILL.md`). In a **parallel** sweep, each
worker starts cold — there's no "previous step" to inherit, because workers don't share state.
Skipping the relaxation step under `distribute` silently produces a different (and usually wrong)
result than the serial version, not just a slower one.

Fix: inside the per-worker function, explicitly relax before measuring:

```python
junction.clearLog()
junction.runSimulation(3e-9, int_step, int_step)   # relax 1-5ns, no logging needed
# ...then run the actual measurement
```

If results still look off after adding relaxation, increase the relax time or decrease `int_step`
before suspecting the physics — both are cheap experiments in a parallel run since the cost is
amortized across workers anyway.
