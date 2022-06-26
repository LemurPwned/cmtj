## `ParallelStack`

### `__init__(self, junctionList: List[cmtj.Junction])`

Initialises a parallel connection of junctions.

### `clearLogs(self)`

Clear all the logs, both of the stack and the junctions
that constitute the stack.

### `getLog(self, junctionId: int)`

Get the logs of a specific junction -- integer id
from the `junctionList`.

### `getLog(self)`

Get the logs of the stack

### `runSimulation(self, totalTime: float, timeStep: float = ..., writeFrequency: float = ...)`

Run the simulation of the stack.

#### **Parameters**

| Name            | Type    | Description                                                                      | Default |
| --------------- | ------- | -------------------------------------------------------------------------------- | ------- |
| **`totalTime`** | `float` | total time of a simulation, give it in seconds. Typical length is in ~couple ns. | `-`     |
| **`timeStep`**  | `float` | the integration step of the RK45 method. Default is 1e-13                        | `...`   |

### `setCoupledCurrentDriver(self, driver: cmtj.ScalarDriver)`

Sets a global current driver for all junctions inside the stack.
Keep in mind the current passed down the stack will be modified
by the coupling constant.

### `setCouplingStrength(self, coupling: float)`

Coupling constant that represents the energy losses as the current
passes through the stack.

### `setExternalFieldDriver(self, driver: cmtj.AxialDriver)`

Sets a external field current driver for all junctions inside the stack.

### `setMagnetistation(self, juncionId: int, layerId: str, mag: cmtj.CVector)`

Set magnetisation on a specific layer in a specific junction.

#### **Parameters**

| Name             | Type  | Description                                         | Default |
| ---------------- | ----- | --------------------------------------------------- | ------- |
| **`junctionId`** | `-`   | the id of the junction (int) as passed in the init. | `-`     |
| **`layerId`**    | `str` | the string id of the layer in the junction.         | `-`     |

## `SeriesStack`

### `__init__(self, junctionList: List[cmtj.Junction])`

Initialises a series connection of junctions.

### `clearLogs(self)`

Clear all the logs, both of the stack and the junctions
that constitute the stack.

### `getLog(self, junctionId: int)`

Get the logs of a specific junction -- integer id
from the `junctionList`.

### `getLog(self)`

Get the logs of the stack

### `runSimulation(self, totalTime: float, timeStep: float = ..., writeFrequency: float = ...)`

Run the simulation of the stack.

#### **Parameters**

| Name            | Type    | Description                                                                      | Default |
| --------------- | ------- | -------------------------------------------------------------------------------- | ------- |
| **`totalTime`** | `float` | total time of a simulation, give it in seconds. Typical length is in ~couple ns. | `-`     |
| **`timeStep`**  | `float` | the integration step of the RK45 method. Default is 1e-13                        | `...`   |

### `setCoupledCurrentDriver(self, driver: cmtj.ScalarDriver)`

Sets a global current driver for all junctions inside the stack.
Keep in mind the current passed down the stack will be modified
by the coupling constant.

### `setCouplingStrength(self, coupling: float)`

Coupling constant that represents the energy losses as the current
passes through the stack.

### `setExternalFieldDriver(self, driver: cmtj.AxialDriver)`

Sets a external field current driver for all junctions inside the stack.

### `setMagnetistation(self, juncionId: int, layerId: str, mag: cmtj.CVector)`

Set magnetisation on a specific layer in a specific junction.

#### **Parameters**

| Name             | Type  | Description                                         | Default |
| ---------------- | ----- | --------------------------------------------------- | ------- |
| **`junctionId`** | `-`   | the id of the junction (int) as passed in the init. | `-`     |
| **`layerId`**    | `str` | the string id of the layer in the junction.         | `-`     |
