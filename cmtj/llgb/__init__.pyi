from typing import Dict, List, Tuple

import cmtj

class LLGBJunction:
    def __init__(self, layers: List[LLGBLayer]) -> None: ...
    def clearLog(self) -> None: ...
    def getLog(self) -> Dict[str, List[float]]: ...
    def runSimulation(
        self,
        totalTime: float,
        timeStep: float = ...,
        writeFrequency: float = ...,
        log: bool = ...,
        solverMode: cmtj.SolverMode = ...,
    ) -> None: ...
    def saveLogs(self, arg0: str) -> None: ...
    def setLayerExternalFieldDriver(
        self, arg0: str, arg1: cmtj.AxialDriver
    ) -> None: ...
    def setLayerTemperatureDriver(self, arg0: str, arg1: cmtj.ScalarDriver) -> None: ...

class LLGBLayer:
    def __init__(
        self,
        id: str,
        mag: cmtj.CVector,
        anis: cmtj.CVector,
        Ms: float,
        thickness: float,
        cellSurface: float,
        demagTensor: List[cmtj.CVector],
        damping: float,
        Tc: float,
        susceptibility: float,
        me: float,
    ) -> None: ...
    def setAnisotropyDriver(self, arg0: cmtj.ScalarDriver) -> None: ...
    def setExternalFieldDriver(self, arg0: cmtj.AxialDriver) -> None: ...
    def setTemperatureDriver(self, arg0: cmtj.ScalarDriver) -> None: ...

def MFAWeissCurie(
    me: float,
    T: float,
    J0: float,
    relax: float = ...,
    tolerance: float = ...,
    maxIter: int = ...,
) -> Tuple[float, float]: ...
def langevin(arg0: float) -> float: ...
def langevinDerivative(arg0: float) -> float: ...
