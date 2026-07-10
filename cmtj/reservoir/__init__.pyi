from typing import Callable, overload

import cmtj

class GroupInteraction:
    def __init__(
        self,
        coordinateMatrix: list[cmtj.CVector],
        junctionList: list[cmtj.Junction],
        topId: str = "free",
    ) -> None:
        """Initialize GroupInteraction for coupled junctions.

        :param coordinateMatrix: List of position vectors for each junction
        :param junctionList: List of junctions to couple
        :param topId: ID of the top layer to use for interactions (default: "free")
        :raises RuntimeError: If coordinate and junction lists have different sizes or are empty
        """
        ...

    def clearLogs(self) -> None:
        """Clear the logs"""
        ...

    @overload
    def getLog(self, junctionIndex: int) -> dict[str, list[float]]:
        """Get the logs for a specific junction.

        :param junctionIndex: Index of the junction
        :raises RuntimeError: If junction index is out of bounds
        :return: Dictionary containing log data
        """
        ...

    @overload
    def getLog(self) -> dict[str, list[float]]:
        """Get the logs for all junctions.

        :return: Dictionary containing log data
        """
        ...

    def runSimulation(self, totalTime: float, timeStep: float = 1e-13, writeFrequency: float = 1e-13) -> None:
        """Run the coupled simulation.

        :param totalTime: Total simulation time
        :param timeStep: Time step for integration
        :param writeFrequency: How often to write data to logs
        :raises RuntimeError: If timeStep > writeFrequency or junctions have incompatible solver modes
        """
        ...

    def setInteractionFunction(
        self,
        function: Callable[[cmtj.CVector, cmtj.CVector, cmtj.Layer, cmtj.Layer], cmtj.CVector],
    ) -> None:
        """Set the interaction function for the coupled junctions.

        :param function: Interaction function.
            Either `computeDipoleInteraction` or `computeDipoleInteractionNoumra` or `nullDipoleInteraction`
            or provide your own custom function.
        """
        ...

def computeDipoleInteraction(
    r1: cmtj.CVector, r2: cmtj.CVector, layer1: cmtj.Layer, layer2: cmtj.Layer
) -> cmtj.CVector:
    """Compute dipole interaction between two junctions (Kanao et al. 2019 PRA).

    :param r1: Position vector of the first junction
    :param r2: Position vector of the second junction
    :param layer1: Magnetic layer of the first junction
    :param layer2: Magnetic layer of the second junction

    :return: Dipole interaction vector
    """
    ...

def computeDipoleInteractionNoumra(
    r1: cmtj.CVector, r2: cmtj.CVector, layer1: cmtj.Layer, layer2: cmtj.Layer
) -> cmtj.CVector:
    """Compute dipole interaction between two junctions (Nomura et al. 2019 JJAP).

    :param r1: Position vector of the first junction
    :param r2: Position vector of the second junction
    :param layer1: Magnetic layer of the first junction
    :param layer2: Magnetic layer of the second junction

    :return: Dipole interaction vector
    """
    ...

def nullDipoleInteraction(r1: cmtj.CVector, r2: cmtj.CVector, layer1: cmtj.Layer, layer2: cmtj.Layer) -> cmtj.CVector:
    """Compute null dipole interaction between two junctions.
    This is a placeholder function that returns a zero vector.

    :param r1: Position vector of the first junction
    :param r2: Position vector of the second junction
    :param layer1: Magnetic layer of the first junction
    :param layer2: Magnetic layer of the second junction

    :return: Zero vector
    """
    ...
