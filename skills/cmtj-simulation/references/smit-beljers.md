# Smit-Beljers static model — `cmtj.models.general_sb`

Full derivation: `docs/physics/sb-model.md`. Use this model only when statics (equilibrium M(H)/R(H)
and FMR frequency without full dynamics) are explicitly wanted — default to the dynamic LLGS
`Junction` API otherwise (see main `SKILL.md`).

## The one unit exception in the whole codebase

`Ms` here is **A/m**, not Tesla, and coordinates are **spherical** (`theta`, `phi` in radians),
not Cartesian. This is the single inversion of the rule everywhere else in CMTJ. Converting from a
Tesla value used elsewhere in the same script:

```python
from cmtj.utils import mu0
Ms1 = 1.0 / mu0    # 1 T -> A/m for the SB model
```

## Minimal pattern

```python
import numpy as np
from collections import defaultdict
from cmtj.models import LayerSB, VectorObj, Solver
from cmtj.utils import mu0

layerA = LayerSB(_id=0, thickness=1e-9, Kv=VectorObj(np.deg2rad(0.), np.deg2rad(0), 1e1), Ks=3e4, Ms=1.0 / mu0)
layerB = LayerSB(_id=1, thickness=1.3e-9, Kv=VectorObj(np.deg2rad(0.), np.deg2rad(0), 1e4), Ks=1e1, Ms=1.2 / mu0)

current_position = [np.deg2rad(89), np.deg2rad(0.1), np.deg2rad(180), np.deg2rad(0.1)]  # theta1,phi1,theta2,phi2 guess
result = defaultdict(list)
for Hmag in Hspace:
    solver = Solver(layers=[layerA, layerB], J1=[1e-4], J2=[0.0], H=VectorObj(np.deg2rad(89), np.deg2rad(0.1), Hmag))
    (t1, p1, t2, p2), frequencies = solver.solve(init_position=current_position)  # frequencies already in GHz
    current_position = [t1, p1, t2, p2]   # reuse as next step's initial guess -- same speed/convergence trick as the dynamic scans
```

Notes:

- `Kv` only uses the `phi` angle of the `VectorObj` (in-plane anisotropy direction) — `theta` is ignored for that field, don't rely on it.
- `frequencies` come back already converted to GHz; don't re-scale them.
- The equilibrium solve is gradient descent (Adam) + a greedy root search on the Hessian determinant — fast for GHz/MHz-scale problems but it's a local solver: a bad `init_position` guess can converge to the wrong minimum, especially near a switching field. Align the initial guess with the external field direction when unsure, and always carry the previous step's result forward during a sweep (as above) rather than resetting the guess every iteration.
- `Solver(..., prefer_numerical_roots=False)` can speed up the root search — see `docs/tipsandtricks.md`.
