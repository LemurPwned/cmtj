---
author:
  - LemurPwned
date: July 2024
title: Tips and tricks
---


# Drivers API

The drivers API allows you to freely control the excitation types in `cmtj` library.
There are two types of drivers:

- `ScalarDriver` -- is a function over time that returns a scalar value.
- `AxialDriver` -- is a driver that returns a 3-vector value. It is a composition of three scalar drivers, one for each component of the vector.

The library provides a few built-in drivers that can be used to define the scalar values. The built-in `ScalarDrivers` drivers that can also be easily used to define the `AxialDrivers`. The built-in drivers are:

- `constantDriver`: A driver that returns a constant value at each time step.
- `sineDriver`: A driver that returns a sinusoidal value at each time step.
- `gaussianImpulseDriver`: A driver that returns a Gaussian impulse at a given time.
- `posSineDriver`: A driver that returns a positive sinusoidal value at each time step.
- `pulseDriver`: A driver that returns a pulse at a given time.
- `stepDriver`: A driver that returns a step function at a given time.
- `trapezoidDriver`: A driver that returns a trapezoidal function at a given time.
- `NullDriver`: A driver that returns zero at each time step (no-op driver)

For more details on the driver parameters see the binding file [here](https://github.com/LemurPwned/cmtj/blob/master/cmtj/__init__.pyi#L14) or [the documentation](../core/#cmtj.constantDriver).
<!-- http://127.0.0.1:8000/cmtj/api/drivers/core/#cmtj.constantDriver
http://127.0.0.1:8000/cmtj/api/core/ -->
## How to define your own drivers?

You can define your own drivers by inheriting from the `ScalarDriver` class. This class has a single method `getCurrentScalarValue` which you need to implement. This method should return the scalar value of the driver at the given time. The time is given in seconds. The driver can be used in the same way as the built-in drivers. Here is an example of a driver that returns a random value at each time step:

```python
from cmtj import (
    ScalarDriver,
    Layer,
    Junction,
    CVector,
    constantDriver,
    AxialDriver,
    NullDriver,
)
import numpy as np


def my_custom_function(time: float) -> float:
    return time * np.random.choice([-1, 1])

# Create a driver with this function
driver = ScalarDriver.getCustomDriver(my_custom_function)

demag = [CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 1)]
layer = Layer(
    "free",
    mag=CVector(0.1, 0.1, 0.9),
    anis=CVector(0.0, 0.0, 1.0),
    Ms=1.0,
    thickness=3e-9,
    cellSurface=0,
    demagTensor=demag,
    damping=3e-3,
)
layer.setReferenceLayer(CVector(0, 0, 1))
junction = Junction([layer], 100, 200)
junction.setLayerExternalFieldDriver(
    "all", AxialDriver(driver, NullDriver(), NullDriver())
)
junction.setLayerAnisotropyDriver("all", constantDriver(150e3))
junction.runSimulation(30e-9, 1e-13, 1e-13)

```

After you've defined the driver these work just like any other driver.
