## `LLGBJunction`
### `__init__(self, layers: List[LLGBLayer])`






### `clearLog(self)`






### `getLog(self)`






### `runSimulation(self,totalTime: float,timeStep: float = ...,writeFrequency: float = ...,log: bool = ...,solverMode: cmtj.SolverMode = ...,)`






### `saveLogs(self, arg0: str)`






### `setLayerExternalFieldDriver(self, arg0: str, arg1: cmtj.AxialDriver)`






### `setLayerTemperatureDriver(self, arg0: str, arg1: cmtj.ScalarDriver)`





  
## `LLGBLayer`
### `__init__(self,id: str,mag: cmtj.CVector,anis: cmtj.CVector,Ms: float,thickness: float,cellSurface: float,demagTensor: List[cmtj.CVector],damping: float,Tc: float,susceptibility: float,me: float,)`






### `setAnisotropyDriver(self, arg0: cmtj.ScalarDriver)`






### `setExternalFieldDriver(self, arg0: cmtj.AxialDriver)`






### `setTemperatureDriver(self, arg0: cmtj.ScalarDriver)`






### `MFAWeissCurie(me: float,T: float,J0: float,relax: float = ...,tolerance: float = ...,maxIter: int = ...,)`






### `langevin(arg0: float)`






### `langevinDerivative(arg0: float)`





  
