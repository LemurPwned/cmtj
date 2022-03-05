## `AxialDriver`

## `Axis`

## `CVector`

## `Junction`

### `__init__(self, layers: List[Layer], filename: str = ...) -> None:...@overloaddef __init__(self, layers: List[Layer], filename: str, Rp: float = ..., Rap: float = ...) -> None:...@overloaddef __init__(self, layers: List[Layer], filename: str, Rx0: List[float], Ry0: List[float], AMR_X: List[float], AMR_Y: List[float], SMR_X: List[float], SMR_Y: List[float], AHE: List[float])`

Creates a junction with a STRIP magnetoresistance.
Each of the Rx0, Ry, AMR, AMR and SMR is list matching the
length of the layers passed (they directly correspond to each layer).
Calculates the magnetoresistance as per: **see reference**:
Spin Hall magnetoresistance in metallic bilayers by Kim, J. et al.

### `__init__(*args, **kwargs) -> Any:...def clearLog(self)`

Reset current simulation state`

### `getLayerMagnetisation(self, layer_id: str) -> CVector:...def getLog(self)`

Retrieve the simulation log [data].

### `getMagnetoresistance(self) -> List[float]:...def runSimulation(self, totalTime: float, timeStep: float = ..., writeFrequency: float = ..., persist: bool = ..., log: bool = ..., calculateEnergies: bool = ...)`

Main run simulation function. Use it to run the simulation.

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

## `Layer`

### `__init__(self, id: str, mag: CVector, anis: CVector, Ms: float, thickness: float, cellSurface: float, demagTensor: List[CVector], dipoleTensor: List[CVector], temperature: float = ..., damping: float = ...)`

The basic structure is a magnetic layer.
Its parameters are defined by the constructor and may be altered
by the drivers during the simulation time.
If you want STT, remember to set the reference vector for the polarisation of the layer.
Use `setReferenceLayer` function to do that.

#### **Parameters**

| Name               | Type            | Description                                                                          | Default |
| ------------------ | --------------- | ------------------------------------------------------------------------------------ | ------- |
| **`id`**           | `str`           | identifiable name for a layer -- e.g. "bottom" or "free".                            | `-`     |
| **`mag`**          | `CVector`       | initial magnetisation. Must be normalised (norm of 1). Used for quicker convergence. | `-`     |
| **`anis`**         | `CVector`       | anisotropy of the layer. A normalised vector                                         | `-`     |
| **`Ms`**           | `float`         | magnetisation saturation. Unit: Tesla [T].                                           | `-`     |
| **`thickness`**    | `float`         | thickness of the layer. Unit: meter [m].                                             | `-`     |
| **`cellSurface`**  | `float`         | surface of the layer, for volume calculation. Unit: meter^2 [m^2].                   | `-`     |
| **`demagTensor`**  | `List[CVector]` | demagnetisation tensor of the layer.                                                 | `-`     |
| **`dipoleTensor`** | `List[CVector]` | dipole tensor of the layer.                                                          | `-`     |

### `createSOTLayer(id: str, mag: CVector, anis: CVector, Ms: float, thickness: float, cellSurface: float, demagTensor: List[CVector], dipoleTensor: List[CVector], damping: float = 0.11, fieldLikeTorque: float = 0, dampingLikeTorque: float = 0)`

Create SOT layer -- including damping and field-like torques that are
calculated based on the effective Spin Hall angles.

#### **Parameters**

| Name               | Type            | Description                                                                          | Default |
| ------------------ | --------------- | ------------------------------------------------------------------------------------ | ------- |
| **`id`**           | `str`           | identifiable name for a layer -- e.g. "bottom" or "free".                            | `-`     |
| **`mag`**          | `CVector`       | initial magnetisation. Must be normalised (norm of 1). Used for quicker convergence. | `-`     |
| **`anis`**         | `CVector`       | anisotropy of the layer. A normalised vector                                         | `-`     |
| **`Ms`**           | `float`         | magnetisation saturation. Unit: Tesla [T].                                           | `-`     |
| **`thickness`**    | `float`         | thickness of the layer. Unit: meter [m].                                             | `-`     |
| **`cellSurface`**  | `float`         | surface of the layer, for volume calculation. Unit: meter^2 [m^2].                   | `-`     |
| **`demagTensor`**  | `List[CVector]` | demagnetisation tensor of the layer.                                                 | `-`     |
| **`dipoleTensor`** | `List[CVector]` | dipole tensor of the layer.                                                          | `-`     |
| **`temperature`**  | `-`             | resting temperature of the layer. Unit: Kelvin [K].                                  | `-`     |

### `createSTTLayer(id: str, mag: CVector, anis: CVector, Ms: float, thickness: float, cellSurface: float, demagTensor: List[CVector], dipoleTensor: List[CVector], damping: float = 0.011, SlonczewskiSpacerLayerParameter: float = 1.0, beta: float = 0.0, spinPolarisation: float = 0.0)`

Create STT layer -- with the standard Slomczewski formulation.

#### **Parameters**

| Name                                  | Type            | Description                                                                                    | Default |
| ------------------------------------- | --------------- | ---------------------------------------------------------------------------------------------- | ------- |
| **`id`**                              | `str`           | identifiable name for a layer -- e.g. "bottom" or "free".                                      | `-`     |
| **`mag`**                             | `CVector`       | initial magnetisation. Must be normalised (norm of 1). Used for quicker convergence.           | `-`     |
| **`anis`**                            | `CVector`       | anisotropy of the layer. A normalised vector                                                   | `-`     |
| **`Ms`**                              | `float`         | magnetisation saturation. Unit: Tesla [T].                                                     | `-`     |
| **`thickness`**                       | `float`         | thickness of the layer. Unit: meter [m].                                                       | `-`     |
| **`cellSurface`**                     | `float`         | surface of the layer, for volume calculation. Unit: meter^2 [m^2].                             | `-`     |
| **`demagTensor`**                     | `List[CVector]` | demagnetisation tensor of the layer.                                                           | `-`     |
| **`dipoleTensor`**                    | `List[CVector]` | dipole tensor of the layer.                                                                    | `-`     |
| **`damping`**                         | `float`         | often marked as alpha in the LLG equation. Damping of the layer. Default 0.011. Dimensionless. | `0.011` |
| **`SlonczewskiSpacerLayerParameter`** | `float`         | Slomczewski parameter. Often marked as lambda.                                                 | `1.0`   |
| **`beta`**                            | `float`         | beta parameter that scales FL/DL ratio.                                                        | `0.0`   |

## `NullDriver(ScalarDriver)`

### `__init__(self)`

An empty driver that does nothing. Use in Axial Driver when
the axis is to be id.

## `ScalarDriver`

### `__init__(self, *args, **kwargs) -> None:...@staticmethoddef getConstantDriver(constantValue: float)`

Constant driver produces a constant signal of a fixed amplitude.

### `getPulseDriver(constantValue: float, amplitude: 'ScalarDriver', period: float, cycle: float)`

Produces a square pulse of certain period and cycle

#### **Parameters**

| Name                | Type             | Description                                                                | Default |
| ------------------- | ---------------- | -------------------------------------------------------------------------- | ------- |
| **`constantValue`** | `float`          | offset (vertical) of the pulse. The pulse amplitude will be added to this. | `-`     |
| **`amplitude`**     | `'ScalarDriver'` | amplitude of the pulse signal                                              | `-`     |
| **`period`**        | `float`          | period of the signal in seconds                                            | `-`     |

### `getSineDriver(constantValue: float, amplitude: 'ScalarDriver', frequency: float, phase: float)`

Produces a sinusoidal signal with some offset (constantValue), amplitude frequency and phase offset.

#### **Parameters**

| Name                | Type             | Description                                                 | Default |
| ------------------- | ---------------- | ----------------------------------------------------------- | ------- |
| **`constantValue`** | `float`          | vertical offset. The sine will oscillate around this value. | `-`     |
| **`amplitude`**     | `'ScalarDriver'` | amplitude of the sine wave                                  | `-`     |
| **`frequency`**     | `float`          | frequency of the sine                                       | `-`     |

### `getStepDriver(constantValue: float, amplitude: float, timeStart: float, timeStop: float)`

Get a step driver. It has amplitude between timeStart and timeStop and 0 elsewhere

#### **Parameters**

| Name                | Type    | Description                                         | Default |
| ------------------- | ------- | --------------------------------------------------- | ------- |
| **`constantValue`** | `float` | offset of the pulse (vertical)                      | `-`     |
| **`amplitude`**     | `float` | amplitude that is added on top of the constantValue | `-`     |
| **`timeStart`**     | `float` | start of the pulse                                  | `-`     |
