## `Axis`

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

## `AxialDriver`

Requires three scalar drivers, one for each axis.

### `__init__(self, x: 'ScalarDriver', y: 'ScalarDriver', z: 'ScalarDriver')`
