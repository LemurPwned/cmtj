## `AxialDriver`

## `Axis`

## `CVector`

## `Junction`

### `__init__(self, layers: List[Layer], filename: str = ...) -> None: ...@overloaddef __init__(self, layers: List[Layer], filename: str, Rp: float = ..., Rap: float = ...) -> None: ...@overloaddef __init__(self,layers: List[Layer],filename: str,Rx0: List[float],Ry0: List[float],AMR_X: List[float],AMR_Y: List[float],SMR_X: List[float],SMR_Y: List[float],AHE: List[float],)`

Creates a junction with a STRIP magnetoresistance.
Each of the Rx0, Ry, AMR, AMR and SMR is list matching the
length of the layers passed (they directly correspond to each layer).
Calculates the magnetoresistance as per: **see reference**:
Spin Hall magnetoresistance in metallic bilayers by Kim, J. et al.

### `__init__(*args, **kwargs) -> Any: ...def clearLog(self)`

Reset current simulation state`

### `getLayerMagnetisation(self, layer_id: str) -> CVector: ...def getLog(self)`

Retrieve the simulation log [data].

### `getMagnetoresistance(self) -> List[float]: ...def runSimulation(self,totalTime: float,timeStep: float = ...,writeFrequency: float = ...,persist: bool = ...,log: bool = ...,calculateEnergies: bool = ...,)`

Main run simulation function.
Use it to run the simulation.

#### **Parameters**

| Name                 | Type | Description                                                                            | Default |
| -------------------- | ---- | -------------------------------------------------------------------------------------- | ------- |
| **`totalTime`**      | `-`  | total time of a simulation, give it in seconds. Typical length is in ~couple ns.       | `-`     |
| **`timeStep`**       | `-`  | the integration step of the RK45 method. Default is 1e-13                              | `-`     |
| **`writeFrequency`** | `-`  | how often is the log saved to? Must be no smaller than `timeStep`. Default is 1e-11.   | `-`     |
| **`persist`**        | `-`  | whether to save to the filename specified in the Junction constructor. Default is true | `-`     |
| **`log`**            | `-`  | if you want some verbosity like timing the simulation. Default is false                | `-`     |

### `setIECDriver(self, bottom_layer: str, top_layer: str, driver: ScalarDriver)`

Set IEC interaction between two layers.
The names of the params are only for convention. The IEC will be set
between bottomLyaer or topLayer, order is irrelevant.

#### **Parameters**

| Name              | Type | Description        | Default |
| ----------------- | ---- | ------------------ | ------- |
| **`bottomLayer`** | `-`  | the first layer id | `-`     |

### `setQuadIECDriver(self, bottom_layer: str, top_layer: str, driver: ScalarDriver)`

Set secondary (biquadratic term) IEC interaction between two layers.
The names of the params are only for convention. The IEC will be set
between bottomLyaer or topLayer, order is irrelevant.

#### **Parameters**

| Name              | Type | Description        | Default |
| ----------------- | ---- | ------------------ | ------- |
| **`bottomLayer`** | `-`  | the first layer id | `-`     |

### `setLayerTemperatureDriver(self, layer_id: str, driver: ScalarDriver) -> None: ...def setLayerAnisotropyDriver(self, layer_id: str, driver: ScalarDriver) -> None: ...def setLayerCurrentDriver(self, layer_id: str, driver: ScalarDriver) -> None: ...def setLayerExternalFieldDriver(self, layer_id: str, driver: AxialDriver) -> None: ...def setLayerMagnetisation(self, layer_id: str, mag: CVector) -> None: ...@overloaddef setLayerOerstedFieldDriver(self, layer_id: str, driver: AxialDriver) -> None: ...def setLayerDampingLikeTorqueDriver(self, layer_id: str, driver: ScalarDriver)`

Set the damping like torque driver for a layer.

#### **Parameters**

| Name           | Type  | Description  | Default |
| -------------- | ----- | ------------ | ------- |
| **`layer_id`** | `str` | the layer id | `-`     |

### `setLayerFieldLikeTorqueDriver(self, layer_id: str, driver: ScalarDriver)`

Set the field like torque driver for a layer.

#### **Parameters**

| Name           | Type  | Description  | Default |
| -------------- | ----- | ------------ | ------- |
| **`layer_id`** | `str` | the layer id | `-`     |

### `setLayerOneFNoise(self, layer_id: str, sources: int, bias: float, scale: float)`

Set 1/f noise for a layer.

#### **Parameters**

| Name           | Type    | Description                                                           | Default |
| -------------- | ------- | --------------------------------------------------------------------- | ------- |
| **`layer_id`** | `str`   | the layer id                                                          | `-`     |
| **`sources`**  | `int`   | the number of generation sources (the more the slower, but more acc.) | `-`     |
| **`bias`**     | `float` | the bias of the noise (p in the Multinomial distribution)             | `-`     |

## `Layer`

### `__init__(self,id: str,mag: CVector,anis: CVector,Ms: float,thickness: float,cellSurface: float,demagTensor: List[CVector],temperature: float = ...,damping: float = ...,)`

The basic structure is a magnetic layer.
Its parameters are defined by the constructor and may be altered
by the drivers during the simulation time.
If you want STT, remember to set the reference vector for the polarisation of the layer.
Use `setReferenceLayer` function to do that.

#### **Parameters**

| Name               | Type      | Description                                                                          | Default |
| ------------------ | --------- | ------------------------------------------------------------------------------------ | ------- |
| **`id`**           | `str`     | identifiable name for a layer -- e.g. "bottom" or "free".                            | `-`     |
| **`mag`**          | `CVector` | initial magnetisation. Must be normalised (norm of 1). Used for quicker convergence. | `-`     |
| **`anis`**         | `CVector` | anisotropy of the layer. A normalised vector                                         | `-`     |
| **`Ms`**           | `float`   | magnetisation saturation. Unit: Tesla [T].                                           | `-`     |
| **`thickness`**    | `float`   | thickness of the layer. Unit: meter [m].                                             | `-`     |
| **`cellSurface`**  | `float`   | surface of the layer, for volume calculation. Unit: meter^2 [m^2].                   | `-`     |
| **`dipoleTensor`** | `-`       | dipole tensor of the layer.                                                          | `-`     |

### `createSOTLayer(id: str,mag: CVector,anis: CVector,Ms: float,thickness: float,cellSurface: float,demagTensor: List[CVector],damping: float = 0.11,fieldLikeTorque: float = 0,dampingLikeTorque: float = 0,) -> "Layer":"""Create SOT layer -- including damping and field-like torques that arecalculated based on the effective Spin Hall angles.:param id: identifiable name for a layer -- e.g. "bottom" or "free".:param mag: initial magnetisation. Must be normalised (norm of 1). Used for quicker convergence.:param anis: anisotropy of the layer. A normalised vector:param Ms: magnetisation saturation. Unit: Tesla [T].:param thickness: thickness of the layer. Unit: meter [m].:param cellSurface: surface of the layer, for volume calculation. Unit: meter^2 [m^2].:param dipoleTensor: dipole tensor of the layer.:param temperature: resting temperature of the layer. Unit: Kelvin [K].:param damping: often marked as alpha in the LLG equation. Damping of the layer. Default 0.011. Dimensionless."""...@staticmethoddef createSTTLayer(id: str,mag: CVector,anis: CVector,Ms: float,thickness: float,cellSurface: float,demagTensor: List[CVector],damping: float = 0.011,SlonczewskiSpacerLayerParameter: float = 1.0,beta: float = 0.0,spinPolarisation: float = 0.0,) -> "Layer":"""Create STT layer -- with the standard Slomczewski formulation.:param id: identifiable name for a layer -- e.g. "bottom" or "free".:param mag: initial magnetisation. Must be normalised (norm of 1). Used for quicker convergence.:param anis: anisotropy of the layer. A normalised vector:param Ms: magnetisation saturation. Unit: Tesla [T].:param thickness: thickness of the layer. Unit: meter [m].:param cellSurface: surface of the layer, for volume calculation. Unit: meter^2 [m^2].:param dipoleTensor: dipole tensor of the layer.:param damping: often marked as alpha in the LLG equation. Damping of the layer. Default 0.011. Dimensionless.:param SlonczewskiSpacerLayerParameter: Slomczewski parameter. Often marked as lambda.:param beta: beta parameter that scales FL/DL ratio.:param spinPolarisation: the spin effectiveness."""...def setAnisotropyDriver(self, driver: ScalarDriver)`

Set anisotropy driver for the layer.
It's scalar. The axis is determined in the layer constructor

### `setTemperatureDriver(self, driver: ScalarDriver)`

Set a driver for the temperature of the layer.
Automatically changes the solver to Euler-Heun.

### `setExternalFieldDriver(self, driver: AxialDriver) -> None: ...def setMagnetisation(self, mag: CVector) -> None: ...def setOerstedFieldDriver(self, driver: AxialDriver) -> None: ...def setDampingLikeTorqueDriver(self, driver: ScalarDriver)`

Set a driver for the damping like torque of the layer.

### `setFieldLikeTorqueDriver(self, driver: ScalarDriver)`

Set a driver for the field like torque of the layer.

### `setReferenceLayer(self, ref: CVector) -> None: ...@overloaddef setReferenceLayer(self, ref: "Reference") -> None: ...def setTopDipoleTensor(self, tensor: List[CVector])`

Set a dipole tensor from the top layer.

### `setBottomDipoleTensor(self, tensor: List[CVector])`

Set a dipole tensor from the bottom layer.

### `getId(self)`

Get Id of the layer

### `setAlternativeSTT(self, setAlternative: bool)`

Switch to an alternative STT forumulation (Taniguchi et al.)
https://iopscience.iop.org/article/10.7567/APEX.11.013005

### `setKappa(self, kappa: float)`

Set the kappa parameter for the layer -- determines SOT mixing
Hdl \* kappa + Hfl
Allows you to turn off Hdl. Turning Hfl is via beta parameter.

## `NullDriver(ScalarDriver)`

### `__init__(self)`

An empty driver that does nothing. Use in Axial Driver when
the axis is to be id.

## `ScalarDriver`

### `__init__(self, *args, **kwargs) -> None: ...@staticmethoddef getConstantDriver(constantValue: float) -> "ScalarDriver":"""Constant driver produces a constant signal of a fixed amplitude.:param constantValue: constant value of the driver (constant offset/amplitude)"""...@staticmethoddef getPulseDriver(constantValue: float, amplitude: "ScalarDriver", period: float, cycle: float)`

Produces a square pulse of certain period and cycle

#### **Parameters**

| Name                | Type | Description                                                                | Default |
| ------------------- | ---- | -------------------------------------------------------------------------- | ------- |
| **`constantValue`** | `-`  | offset (vertical) of the pulse. The pulse amplitude will be added to this. | `-`     |
| **`amplitude`**     | `-`  | amplitude of the pulse signal                                              | `-`     |
| **`period`**        | `-`  | period of the signal in seconds                                            | `-`     |

### `getSineDriver(constantValue: float, amplitude: "ScalarDriver", frequency: float, phase: float)`

Produces a sinusoidal signal with some offset (constantValue), amplitude frequency and phase offset.

#### **Parameters**

| Name                | Type             | Description                                                 | Default |
| ------------------- | ---------------- | ----------------------------------------------------------- | ------- |
| **`constantValue`** | `float`          | vertical offset. The sine will oscillate around this value. | `-`     |
| **`amplitude`**     | `"ScalarDriver"` | amplitude of the sine wave                                  | `-`     |
| **`frequency`**     | `float`          | frequency of the sine                                       | `-`     |

### `getStepDriver(constantValue: float, amplitude: float, timeStart: float, timeStop: float) -> "ScalarDriver":"""Get a step driver. It has amplitude between timeStart and timeStop and 0 elsewhere:param constantValue: offset of the pulse (vertical):param amplitude: amplitude that is added on top of the constantValue:param timeStart: start of the pulse:param timeStop: when the pulse ends"""...@staticmethoddef getTrapezoidDriver(constantValue: float,amplitude: float,timeStart,edgeTime: float,steadyTime: float,)`

Create Trapezoid driver. Has a rising and a falling edge.

#### **Parameters**

| Name                | Type    | Description                                         | Default |
| ------------------- | ------- | --------------------------------------------------- | ------- |
| **`constantValue`** | `float` | offset of the pulse (vertical)                      | `-`     |
| **`amplitude`**     | `float` | amplitude that is added on top of the constantValue | `-`     |
| **`timeStart`**     | `float` | start of the pulse                                  | `-`     |
| **`edgeTime`**      | `-`     | time it takes to reach the maximum amplitude        | `-`     |

### `getGaussianImpulseDriver(constantValue: float, amplitude: float, t0: float, sigma: float)`

Gaussian impulse driver. It has amplitude starts at t0 and falls off with sigma.

        Formula:
        A * exp(-((t - t0) ** 2) / (2 * sigma ** 2))

#### **Parameters**

| Name                | Type    | Description                                         | Default |
| ------------------- | ------- | --------------------------------------------------- | ------- |
| **`constantValue`** | `float` | offset of the pulse (vertical)                      | `-`     |
| **`amplitude`**     | `float` | amplitude that is added on top of the constantValue | `-`     |
| **`t0`**            | `float` | start of the pulse                                  | `-`     |

### `getGaussianStepDriver(constantValue: float, amplitude: float, t0: float, sigma: float)`

Gaussian step driver (erf function). It has amplitude starts at t0 and falls off with sigma.

        Formula:
        f(t) = constantValue + amplitude * (1 + erf((t - t0) / (sigma * sqrt(2))))

#### **Parameters**

| Name                | Type    | Description                                         | Default |
| ------------------- | ------- | --------------------------------------------------- | ------- |
| **`constantValue`** | `float` | offset of the pulse (vertical)                      | `-`     |
| **`amplitude`**     | `float` | amplitude that is added on top of the constantValue | `-`     |
| **`t0`**            | `float` | start of the pulse                                  | `-`     |

### `getPosSineDriver(constantValue: float, amplitude: float, frequency: float, phase: float)`

Produces a positive sinusoidal signal with some offset (constantValue), amplitude frequency and phase offset.

#### **Parameters**

| Name                | Type    | Description                                                 | Default |
| ------------------- | ------- | ----------------------------------------------------------- | ------- |
| **`constantValue`** | `float` | vertical offset. The sine will oscillate around this value. | `-`     |
| **`amplitude`**     | `float` | amplitude of the sine wave                                  | `-`     |
| **`frequency`**     | `float` | frequency of the sine                                       | `-`     |

### `getPulseDriver(constantValue: float, amplitude: float, period: float, cycle: float)`

Produces a square pulse of certain period and cycle

#### **Parameters**

| Name                | Type    | Description                                                                | Default |
| ------------------- | ------- | -------------------------------------------------------------------------- | ------- |
| **`constantValue`** | `float` | offset (vertical) of the pulse. The pulse amplitude will be added to this. | `-`     |
| **`amplitude`**     | `float` | amplitude of the pulse signal                                              | `-`     |
| **`period`**        | `float` | period of the signal in seconds                                            | `-`     |

## `SolverMode`

## `Reference`
