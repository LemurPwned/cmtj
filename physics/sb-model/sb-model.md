---
author:
  - LemurPwned
date: Feburary 2023
title: Static models -- Smit-Beljers
---
# Smit-Beljers model

## Introduction
The basic model was introduced by Smit and Beljers 1955. However, there were nice extensions of the model, notably, by L. Baselgia, M. Warden, and F. Waldner in the *Derivation of the resonance frequency from the free energy of ferromagnets*, 1988. One of more recent ones by R. L. Rodríguez-Suárez, S. M. Rezende, and A. Azevedo from *Ferromagnetic resonance investigation of the residual coupling in spin-valve systems*, PRB 2005 is the one I prefer to use for explanation.

## The model
The model is capable of simulating the magnetisation in equilibrium and the resonance frequency mode (fmr = ferromagnetic resonance). The former basically allows us to reproduce the hysteresis M(H) or R(H) loops etc. The latter gives use the frequency of the oscillations of the magnetisation in the presence of an external magnetic field. 

For the full model we need:

1. Create energy expression of the system.
2. Obtain equilibirum magnetisation position of the system. 
3. Compute the hessian matrix of the energy expression.
4. Find the roots of the det(hessian). 

Each layer is described by the following parameters:
##### Parameter table
|    Parameter    |                   Description                    |      Units      |
| :-------------: | :----------------------------------------------: | :-------------: |
|    $\theta$     |       Polar angle of magnetisation vector        |       rad       |
|     $\phi$      |     Azimuthal angle of magnetisation vector      |       rad       |
| $t_\textrm{FM}$ |            Thickness of the FM layer             |        m        |
|      $M_s$      |             Magnetisation saturation             |  $\frac{A}{m}$  |
|      $K_s$      | Surface anisotropy (perpendicular, out of plane) | $\frac{J}{m^3}$ |
|      $K_s$      |          Volume anisotropy  (in-plane)           | $\frac{J}{m^3}$ |
|    $\alpha$     |      Azimuthal angle of in-plane anisotropy      | $\frac{J}{m^3}$ |
|      $J_1$      |   Interlayer exchange coupling value (linear)    | $\frac{J}{m^2}$ |
|      $J_2$      |  Interlayer exchange coupling value (quadratic)  | $\frac{J}{m^2}$ |

And the external field is also given by it's own ($\theta$, $\phi$, $H$) where $H$ is the magnitude of the field.

!!! note

    This assumes physical, not mathematical spherical coordinates! Here, polar angle $\theta$ is measured from the positive $z$ axis, and azimuthal angle $\phi$ is measured from the positive $x$ axis.

The energy expression per layer is given by the following:

$$
    \varepsilon_i = - \mu_0 M_s \vec{m}_i \cdot \vec{H} + 
        (- K_s + \frac{1}{2}\mu_0 M_s^2) m_z^2 
        - K_v \vec{m} \cdot \vec{K}_{dir} 
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
The equilibrium magnetisation, $(\theta^*, \phi^*)$, is then used to compute the hessian matrix of the energy expression. The hessian matrix is a matrix of second derivatives of the energy expression. We use it to find the roots of the det(hessian) expression. Those roots designate the frequencies of the resonance mode. 

## Root finding
Root finding algorithm is a native greedy search, but for GHz or MHz frequencies it's pretty fast and precise enough (you can set the tolerance in the parameters).