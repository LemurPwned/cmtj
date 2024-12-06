from typing import overload

import cmtj

class ParallelStack:
    def __init__(
        self,
        junctionList: list[cmtj.Junction],
        topId: str = "free",
        bottomId: str = "bottom",
        phaseOffset: float = 0,
    ) -> None:
        """
        Initialises a parallel connection of junctions.
        Layer ids are used to identify the layers in the junctions and for resistance calculations.
        :param junctionList: list of junctions to be connected in parallel.
        :param topId: the string id of the top layer in the stack. Default is "free".
        :param bottomId: the string id of the bottom layer in the stack. Default is "bottom".
        :param phaseOffset: the phase offset between the junctions. Default is 0.
        """
        ...

    def clearLogs(self) -> None:
        """
        Clear all the logs, both of the stack and the junctions
        that constitute the stack.
        """
        ...

    @overload
    def getLog(self, junctionId: int) -> dict[str, list[float]]:
        """
        Get the logs of a specific junction -- integer id
        from the `junctionList`.
        :param junctionId: integer junction id as was passed in the init.
        """
        ...

    @overload
    def getLog(self) -> dict[str, list[float]]:
        """
        Get the logs of the stack
        """
        ...

    def runSimulation(self, totalTime: float, timeStep: float = ..., writeFrequency: float = ...) -> None:
        """
        Run the simulation of the stack.
        :param totalTime: total time of a simulation, give it in seconds. Typical length is in ~couple ns.
        :param timeStep: the integration step of the RK45 method. Default is 1e-13
        :param writeFrequency: how often is the log saved to? Must be no smaller than `timeStep`. Default is 1e-11.
        """
        ...

    def setCoupledCurrentDriver(self, driver: cmtj.ScalarDriver) -> None:
        """
        Sets a global current driver for all junctions inside the stack.
        Keep in mind the current passed down the stack will be modified
        by the coupling constant.
        :param driver: the current driver to be set.
        """
        ...

    def setCouplingStrength(self, coupling: float) -> None:
        """
        Coupling constant that represents the energy losses as the current
        passes through the stack.
        :param coupling: the coupling strength (or the losses)
        """
        ...

    def setExternalFieldDriver(self, driver: cmtj.AxialDriver) -> None:
        """
        Sets a external field current driver for all junctions inside the stack.
        :param driver: the field driver to be set.
        """
        ...

    def setMagnetisation(self, junctionId: int, layerId: str, mag: cmtj.CVector) -> None:
        """
        Set magnetisation on a specific layer in a specific junction.
        :param junctionId: the id of the junction (int) as passed in the init.
        :param layerId: the string id of the layer in the junction.
        :param mag: the magnetisation to be set.
        """
        ...

    def getMagnetisation(self, junction: int, layerId: str) -> cmtj.CVector:
        """Get the magnetisation of a specific layer in a specific junction.
        :param junction: the id of the junction (int) as passed in the init.
        :param layerId: the string id of the layer in the junction."""
        ...

    def getJunction(self, junctionId: int) -> cmtj.Junction:
        """Get a specific junction from the stack. Returns a reference.
        :param junctionId: the id of the junction (int) as passed in the init.
        """
        ...

class SeriesStack:
    def __init__(
        self,
        junctionList: list[cmtj.Junction],
        topId: str = "free",
        bottomId: str = "bottom",
        phaseOffset: float = 0,
    ) -> None:
        """
        Initialises a series connection of junctions.
        Layer ids are used to identify the layers in the junctions and for resistance calculations.
        :param junctionList: list of junctions to be connected in series.
        :param topId: the string id of the top layer in the stack. Default is "free".
        :param bottomId: the string id of the bottom layer in the stack. Default is "bottom".
        :param phaseOffset: the phase offset between the junctions. Default is 0.
        """
        ...

    def clearLogs(self) -> None:
        """
        Clear all the logs, both of the stack and the junctions
        that constitute the stack.
        """
        ...

    @overload
    def getLog(self, junctionId: int) -> dict[str, list[float]]:
        """
        Get the logs of a specific junction -- integer id
        from the `junctionList`.
        :param junctionId: integer junction id as was passed in the init.
        """
        ...

    @overload
    def getLog(self) -> dict[str, list[float]]:
        """
        Get the logs of the stack
        """
        ...

    def runSimulation(self, totalTime: float, timeStep: float = ..., writeFrequency: float = ...) -> None:
        """
        Run the simulation of the stack.
        :param totalTime: total time of a simulation, give it in seconds. Typical length is in ~couple ns.
        :param timeStep: the integration step of the RK45 method. Default is 1e-13
        :param writeFrequency: how often is the log saved to? Must be no smaller than `timeStep`. Default is 1e-11.
        """
        ...

    def setCoupledCurrentDriver(self, driver: cmtj.ScalarDriver) -> None:
        """
        Sets a global current driver for all junctions inside the stack.
        Keep in mind the current passed down the stack will be modified
        by the coupling constant.
        :param driver: the current driver to be set.
        """
        ...

    @overload
    def setCouplingStrength(self, coupling: float) -> None:
        """
        Coupling constant that represents the energy losses as the current
        passes through the stack.
        :param coupling: the coupling strength (or the losses)
        """
        ...

    @overload
    def setCouplingStrength(self, coupling: list[float]) -> None:
        """
        Coupling constant that represents the energy losses as the current
        passes through the stack.
        :param coupling: the coupling strength (or the losses) for each junction.
            Must be the one less than length of the junction vector, i.e. len(junctionList)-1 .
        """
        ...

    def setExternalFieldDriver(self, driver: cmtj.AxialDriver) -> None:
        """
        Sets a external field current driver for all junctions inside the stack.
        :param driver: the field driver to be set.
        """
        ...

    def setMagnetisation(self, junctionId: int, layerId: str, mag: cmtj.CVector) -> None:
        """
        Set magnetisation on a specific layer in a specific junction.
        :param junctionId: the id of the junction (int) as passed in the init.
        :param layerId: the string id of the layer in the junction.
        :param mag: the magnetisation to be set.
        """
        ...

    def getMagnetisation(self, junction: int, layerId: str) -> cmtj.CVector:
        """Get the magnetisation of a specific layer in a specific junction.
        :param junction: the id of the junction (int) as passed in the init.
        :param layerId: the string id of the layer in the junction."""
        ...

    def getJunction(self, junctionId: int) -> cmtj.Junction:
        """Get a specific junction from the stack. Returns a reference.
        :param junctionId: the id of the junction (int) as passed in the init.
        """
        ...
