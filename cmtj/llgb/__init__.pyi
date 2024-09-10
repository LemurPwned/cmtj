import cmtj

class LLGBJunction:
    """LLGB Junction class."""

    def __init__(self, layers: list[LLGBLayer]) -> None:
        """Initialises a LLGB junction with layers.
        :param layers: list of LLGB layers."""
        ...

    def clearLog(self) -> None:
        """Clears the simulation log of the junction."""
        ...

    def getLog(self) -> dict[str, list[float]]:
        """Returns the simulation log of the junction."""
        ...

    def runSimulation(
        self,
        totalTime: float,
        timeStep: float = ...,
        writeFrequency: float = ...,
        log: bool = ...,
        solverMode: cmtj.SolverMode = ...,
    ) -> None:
        """Runs the simulation of the junction.
        :param totalTime: total simulation time.
        :param timeStep: time step.
        :param writeFrequency: frequency of writing to the log.
        :param log: whether to log the simulation.
        :param solverMode: solver mode.
        """
        ...

    def saveLogs(self, arg0: str) -> None:
        """Saves the simulation logs to a file.
        :param arg0: file path."""
        ...

    def setLayerExternalFieldDriver(self, layerId: str, driver: cmtj.AxialDriver) -> None:
        """Set an external field driver for a layer.
        :param layerId: the id of the layer.
        :param driver: the field driver to be set."""
        ...

    def setLayerTemperatureDriver(self, layerId: str, driver: cmtj.ScalarDriver) -> None:
        """Set a temperature driver for a layer.
        :param layerId: the id of the layer.
        :param driver: the temperature driver to be set.
        """
        ...

class LLGBLayer:
    """LLGB Layer class."""

    def __init__(
        self,
        id: str,
        mag: cmtj.CVector,
        anis: cmtj.CVector,
        Ms: float,
        thickness: float,
        cellSurface: float,
        demagTensor: list[cmtj.CVector],
        damping: float,
        Tc: float,
        susceptibility: float,
        me: float,
    ) -> None:
        """Creates a LLGB layer.
        :param id: layer id.
        :param mag: magnetisation.
        :param anis: anisotropy axis.
        :param Ms: saturation magnetisation.
        :param thickness: thickness.
        :param cellSurface: cell surface.
        :param demagTensor: demagnetisation tensor.
        :param damping: damping factor.
        :param Tc: Curie temperature.
        :param susceptibility: susceptibility.
        :param me: equilibrium magnetisation.
        """
        ...

    def setAnisotropyDriver(self, driver: cmtj.ScalarDriver) -> None:
        """Sets an anisotropy driver.
        :param driver: the anisotropy driver to be set."""
        ...

    def setExternalFieldDriver(self, driver: cmtj.AxialDriver) -> None:
        """Sets an external field driver.
        :param driver: the field driver to be set."""
        ...

    def setTemperatureDriver(self, driver: cmtj.ScalarDriver) -> None:
        """Sets a temperature driver.
        :param driver: the temperature driver to be set."""
        ...

def MFAWeissCurie(
    me: float,
    T: float,
    J0: float,
    relax: float = ...,
    tolerance: float = ...,
    maxIter: int = ...,
) -> tuple[float, float]:
    """Mean Field Approximation for Weiss Curie temperature.
    :param me: equilibrium magnetisation.
    :param T: temperature.
    :param J0: exchange coupling.
    :param relax: relaxation factor.
    :param tolerance: tolerance for convergence.
    :param maxIter: maximum number of iterations."""
    ...

def langevin(arg0: float) -> float: ...
def langevinDerivative(arg0: float) -> float: ...
