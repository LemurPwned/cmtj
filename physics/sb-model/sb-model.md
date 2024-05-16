---
author:
  - LemurPwned
date: Feburary 2023
title: Static models -- Smit-Beljers
---

# Smit-Beljers model

## Introduction

The basic model was introduced by Smit and Beljers 1955. However, there were nice extensions of the model, notably, by L. Baselgia, M. Warden, and F. Waldner in the _Derivation of the resonance frequency from the free energy of ferromagnets_, 1988. One of more recent ones by R. L. Rodríguez-Suárez, S. M. Rezende, and A. Azevedo from _Ferromagnetic resonance investigation of the residual coupling in spin-valve systems_, PRB 2005 is the one I prefer to use for explanation.

## The model

The model is capable of simulating the magnetisation in equilibrium and the resonance frequency mode (fmr = ferromagnetic resonance). The former basically allows us to reproduce the hysteresis M(H) or R(H) loops etc. The latter gives use the frequency of the oscillations of the magnetisation in the presence of an external magnetic field.

For the full model we need:

1. Create energy expression of the system.
2. Obtain equilibirum magnetisation position of the system.
3. Compute the hessian matrix $H$ of the energy expression.
4. Find the roots of the $\det H$.

Each layer is described by the following parameters:

##### Parameter table

|    Parameter    |                  Description                   |      Units      |
| :-------------: | :--------------------------------------------: | :-------------: |
|    $\theta$     |      Polar angle of magnetisation vector       |       rad       |
|     $\phi$      |    Azimuthal angle of magnetisation vector     |       rad       |
| $t_\textrm{FM}$ |           Thickness of the FM layer            |        m        |
|      $M_s$      |            Magnetisation saturation            |  $\frac{A}{m}$  |
|      $K_s$      | Shape anisotropy (perpendicular, out of plane) | $\frac{J}{m^3}$ |
|      $K_v$      |          Volume anisotropy (in-plane)          | $\frac{J}{m^3}$ |
|    $\alpha$     |     Azimuthal angle of in-plane anisotropy     |       rad       |
|      $J_1$      |  Interlayer exchange coupling value (linear)   | $\frac{J}{m^2}$ |
|      $J_2$      | Interlayer exchange coupling value (quadratic) | $\frac{J}{m^2}$ |

And the external field is also given by it's own ($\theta$, $\phi$, $H$) where $H$ is the magnitude of the field.

!!! note

    This assumes physical, not mathematical spherical coordinates! Here, polar angle $\theta$ is measured from the positive $z$ axis, and azimuthal angle $\phi$ is measured from the positive $x$ axis.

The energy expression per layer is given by the following:

$$
    \varepsilon_i = - \mu_0 M_s \vec{m}_i \cdot \vec{H} +
        (- K_s + \frac{1}{2}\mu_0 M_s^2) m_z^2
        - K_v (\vec{m} \cdot \vec{K}_{dir})^2
        - J_1 \vec{m}_i \cdot \vec{m}_{i+1}
        - J_2 (\vec{m}_i \cdot \vec{m}_{i+1})^2
$$

where $\vec{m}$ is in the spherical coordinates, and $\vec{K}_{dir}$ is the direction of the in-plane anisotropy.

$$
    \vec{m} = \begin{bmatrix}
        \sin\theta\cos\phi \\
        \sin\theta\sin\phi \\
        \cos\theta          \\
    \end{bmatrix}
$$

Similarly, $H$ for its respective coordinates.
Finally, $\vec{K}_{dir}$ is given by:

$$
\vec{K}_{dir} = \begin{bmatrix}
    \cos\alpha \\
    \sin\alpha \\
    0          \\
\end{bmatrix}
$$

## Solution

## Equilibrium magnetisation

The solution for equilibrium magnetisation is solved using the Adam gradient descent. It is a quick and reliable algorithm for finding energy minima. It operates on the first derivative of the $\varepsilon_i$ expression, i.e. $\frac{\partial\varepsilon}{\partial\theta}$ and $\frac{\partial\varepsilon}{\partial\phi}$. If the gradient is found within a certain tolerance, the algorithm stops. Otherwise, it continues to iterate until the maximum number of iterations is reached.

## Resonance frequency

The equilibrium magnetisation, $(\theta^*, \phi^*)$, is then used to compute the hessian matrix of the energy expression. The hessian matrix is a matrix of second derivatives of the energy expression. We use it to find the roots of the $\det H$(hessian) expression. Those roots designate the frequencies of the resonance mode.

## Root finding

Root finding algorithm is a naive greedy search, but for GHz or MHz frequencies it's pretty fast and precise enough (you can set the tolerance in the parameters).

## Runnning the model

The same example can be found in the [`Examples` section](../tutorials/SBModel.ipynb), expanded with dynamic spherical approach.
Below is an example of how the model can be used, based on a system with 2 ferromagnetic layers:

```python
import numpy as np

from collections import defaultdict
from cmtj.models import LayerSB, VectorObj, Solver
from cmtj.utils import mu0
from tqdm import tqdm

Ms1 = 1. / mu0 # here we pass the saturation magnetisation in A/m, but in the dynamic model we use T!
Ms2 = 1.2 / mu0
layerA = LayerSB(
    _id=0,
    thickness=1e-9,
    Kv=VectorObj(np.deg2rad(0.), np.deg2rad(0), 1e1), # for the Kv only phi angle counts !
    Ks=3e4,
    Ms=Ms1,
)
layerB = LayerSB(
    _id=1,
    thickness=1.3e-9,
    Kv=VectorObj(np.deg2rad(0.), np.deg2rad(0), 1e4),
    Ks=1e1,
    Ms=Ms2,
)

# we indicate the "guess" of the initial position
# it's generally good to align it with the field, but it's not necessary
current_position = [
    np.deg2rad(89),
    np.deg2rad(0.1),
    np.deg2rad(180),
    np.deg2rad(0.1)
]
Hspace = np.linspace(-400e3, 400e3, 100)
result_dictionary = defaultdict(list)
# we perform a sweep over the field magnitude
for Hmag in tqdm(Hspace):
    solver = Solver(
        layers=[layerA, layerB],
        J1=[1e-4],
        J2=[0.],
        H=VectorObj(np.deg2rad(89), np.deg2rad(0.1), Hmag)
    )
    # check for additional parameters in the solver
    # such as gradient convergence tolerance, max iterations, etc.
    # also in the solver there are root finding parameters
    (t1, p1, t2, p2), frequencies = solver.solve(init_position=current_position)
    # frequencies are already in GHz
    for frequency in frequencies:
        result_dictionary["frequency"].append(frequency)
        result_dictionary["Hmag"].append(Hmag)

    # we note the final position of the magnetisation in spherical coordinates
    result_dictionary["theta_1"].append(t1)
    result_dictionary["phi_1"].append(p1)
    result_dictionary["theta_2"].append(t2)
    result_dictionary["phi_2"].append(p2)
    # we reuse the previous solution as the initial guess for the next iteration
    current_position = [t1, p1, t2, p2]
```

## References

1. Rodríguez-Suárez, R. L., Rezende, S. M. & Azevedo, A. Ferromagnetic resonance investigation of the residual coupling in spin-valve systems. Phys. Rev. B 71, 224406 (2005).
2. Baselgia, L. et al. Derivation of the resonance frequency from the free energy of ferromagnets. Phys. Rev. B 38, 2237–2242 (1988).
