---
author:
  - LemurPwned
date: October 2022
title: Parallelism
---

# Parallelism in CMTJ

## Overview

`cmtj` provides a simple way to parallelise simulation sweeps. For instance, while performing the Voltage Spin Diode (VSD) experiments, we need to sweep the frequency and for each frequency value we further need to sweep with the field. For each field and frequency pair we need to compute the VSD DC resistance. Normally, this scales as $O(NK)$ where $N$ is the number of field values and $K$ is the number of frequency values. However, with parallelism, we can reduce this to $O(K)$ if we would be able to process each frequency value in parallel.

In practice, we use worker pools which are approximately equal to the number of threads in the CPU. Each worker is assigned a frequency value and it computes the VSD DC resistance for all the field values. Whenever the worker is finished, it is assigned another frequency value.

## General caveats

!!! note

    Keep in mind that these issues are not specific to CMTJ, but are a general problem of any parallelisation.

Before we move on to a specific example, let's discuss when we should be careful with simulating in parallel.
When we sweep with an external field, which in the experiment is usually continously swept, we need to realise that if we were to compute the experiment for each separate field in a separate process, we would not carry the relaxed state of magnetisation computed in the previous field step to the next. This is of course a natural consequence of putting things in parallel and can be easily remedided by the following precautions:

- relax the magnetisation the magentisation in each parallel process by calling `runSimulation` for a short time, say between 1-5ns (depending on the complexity of the system),
- increase the simulation time, since reaching the stable state may take a little longer,
- decrease the time step (integration step) -- this is more costly computationally, but it is also ensuring absolute convergence.

We generally recommend to use the `runSimulation` function to relax the magnetisation in each parallel process, and then to use the `runSimulation` function to compute the experiment in parallel. This is the most efficient way to do it, since short relaxation period does not incur a large computational cost, and because we paralellise the experiment anyway, that cost is offset by running on multiple threads.

With this in mind, let's move on to a specific example.

## Problem description

Suppose we want to compute the VSD DC resistance for a given frequency and field value, but also check how the VSD spectra differ for different values of the interlayer exchange coupling (IEC) constant `J`.
What we want to effectively do is we want to fun a VSD experiment for each of the `J` values, wherein each experiment is a function of the frequency and field.
Below, we outline the steps to parallelise such a problem.

## Parallelising the VSD experiment

Here's what we need:

- A function that computes the VSD DC resistance for a given frequency and field value.
- A range of frequencies and fields.
- A distributing function that assigns the frequency and field values to the workers.

The last point is taken care of by importing the `distribute` function from `cmtj.utils.parallel`. The `distribute` function takes a function and a list of arguments and distributes the arguments to the workers.

An example function that computes the VSD DC resistance for a given system with a IEC constant `J` and for a specific frequency and field value:

```python
def simulate_vsd(J, H, frequency):
    """
    This function computes the VSD DC resistance for a
    given frequency f and field value H.
    It also takes the current coupling value J as an argument.
    """
    int_step = 5e-13
    demag = [CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 1)]

    Ku = 0.8e3
    alpha = 0.024  # 0.024
    l1 = Layer(
        "free",
        mag=CVector(1, 0.1, 0.1),
        anis=CVector(0, 1., 0.),  # direction of the anisotropy
        Ms=1.03,
        thickness=2.1e-9,
        cellSurface=surf,
        demagTensor=demag,
        damping=alpha)

    l2 = Layer(
        "bottom",
        mag=CVector(1, 0.1, 0.1),
        anis=CVector(1, 0, 0),  # direction of the anisotropy
        Ms=1.65,
        thickness=6e-9,
        cellSurface=surf,
        demagTensor=demag,
        damping=alpha)
    j1 = Junction([l1, l2], 163.5, 176)
    j1.setLayerAnisotropyDriver("free", ScalarDriver.getConstantDriver(Ku))
    j1.setLayerAnisotropyDriver("bottom", ScalarDriver.getConstantDriver(1e12))
    j1.setIECDriver("free", "bottom", ScalarDriver.getConstantDriver(J))
    hangle = FieldScan.angle2vector(90, 90, H)
    j1.clearLog()
    j1.setLayerExternalFieldDriver(
        "all",
        AxialDriver(ScalarDriver.getConstantDriver(hangle.x),
                    ScalarDriver.getConstantDriver(hangle.y),
                    ScalarDriver.getConstantDriver(hangle.z)))
    j1.setLayerOerstedFieldDriver(
        "free",
        AxialDriver(
            NullDriver(),
            ScalarDriver.getSineDriver(0, 5, frequency, 0),
            NullDriver(),
        ))
    j1.runSimulation(60e-9, int_step, int_step, solverMode=SolverMode.RK4)
    log = j1.getLog()
    dynamicR = log['R_free_bottom']

    dynamicI = np.sin(2 * math.pi * frequency * np.asarray(log['time']))
    vmix = compute_sd(dynamicR, dynamicI, int_step)
    return vmix
```

Writing such functions is a good idea in general, since they can be also used in the normal serial processing as well. So a good practice is to write such a function first, test it in serial/iterative sweep and then use `distribute` to parallelise it easily. Here's how we would use the `distribute` function to parallelise the VSD experiment:

```python
from cmtj.utils.parallel import distribute

Js = np.linspace(1e-4, 1e-3, 10)
Hrange = np.linspace(-15e3, 15e3, 100, endpoint=True)
fscan = np.arange(1e9, 6.2e9, 0.2e9)
VSD = np.zeros((len(Js), len(fscan), len(Hrange)), dtype=np.float32)
Js = np.around(Js, decimals=7)
for res in distribute(simulate_vsd, [Js, Hrange, fscan]):
    (k, i, j), out = res
    VSD[k, j, i] = out
```

The `distribute` function returns a generator, which we can iterate over to get the results. The results are returned as a tuple of two elements: the first element is a tuple of the indices of the result, and the second element is the result itself. In our case, the result is the VSD DC resistance, and the indices are the indices of the `J`, `H` and `f` values (in that order specifically, we need to match the function signature of `simulate_vsd` with the second argument to `distribute` function). We can use these indices to assign the result to the correct place in the `VSD` array.

Very simple!
Now we can plot the VSD spectra for different values of the IEC constant `J`:

```python
fig, axs = plt.subplots(len(Js), 1, figsize=(8, 6), dpi=300)
for k in range(len(Js)):
    axs[k].pcolor( Hrange/1e3, fscan / 1e9, VSD[k], shading='auto')
    axs[k].set_xlabel("H (kA/m)")
    axs[k].set_ylabel("Frequency (GHz)")
    axs[k].set_title("J = {}".format(Js[k]))
```
