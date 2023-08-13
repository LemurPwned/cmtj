from typing import Any, ClassVar, Dict, List, overload

xaxis: Axis
yaxis: Axis
zaxis: Axis


def c_dot(arg0: CVector, arg1: CVector) -> float:
    """Compute dot (scalar) product of two CVectors."""
    ...


class AxialDriver:
    @overload
    def __init__(self, x_driver: ScalarDriver, y_driver: ScalarDriver,
                 z_driver: ScalarDriver) -> None:
        ...

    @overload
    def __init__(self, arg0: List[ScalarDriver]) -> None:
        ...

    @overload
    def __init__(*args, **kwargs) -> Any:
        ...

    @overload
    def applyMask(self, arg0: CVector) -> None:
        ...

    @overload
    def applyMask(self, arg0: List[int]) -> None:
        ...

    @overload
    def applyMask(*args, **kwargs) -> Any:
        ...

    def getCurrentAxialDrivers(self, arg0: float) -> CVector:
        ...

    def getVectorAxialDriver(self, arg0: float, arg1: float) -> AxialDriver:
        ...


class Axis:
    __entries: Any = ...
    xaxis: Any = ...
    yaxis: Any = ...
    zaxis: Any = ...

    def __init__(self, value: int) -> None:
        ...

    def __eq__(self, other: object) -> bool:
        ...

    def __getstate__(self) -> int:
        ...

    def __hash__(self) -> int:
        ...

    def __index__(self) -> int:
        ...

    def __int__(self) -> int:
        ...

    def __ne__(self, other: object) -> bool:
        ...

    def __setstate__(self, state: int) -> None:
        ...

    @property
    def name(self) -> Any:
        ...

    @property
    def __doc__(self) -> Any:
        ...

    @property
    def __members__(self) -> Any:
        ...


class CVector:
    def __init__(self, x: float, y: float, z: float) -> None:
        ...

    def length(self) -> float:
        ...

    @property
    def x(self) -> float:
        ...

    @x.setter
    def x(self, val: float) -> None:
        ...

    @property
    def y(self) -> float:
        ...

    @y.setter
    def y(self, val: float) -> None:
        ...

    @property
    def z(self) -> float:
        ...

    @z.setter
    def z(self, val: float) -> None:
        ...


class Junction:
    @overload
    def __init__(self, layers: List[Layer], filename: str = ...) -> None:
        ...

    @overload
    def __init__(self,
                 layers: List[Layer],
                 filename: str,
                 Rp: float = ...,
                 Rap: float = ...) -> None:
        ...

    @overload
    def __init__(
        self,
        layers: List[Layer],
        filename: str,
        Rx0: List[float],
        Ry0: List[float],
        AMR_X: List[float],
        AMR_Y: List[float],
        SMR_X: List[float],
        SMR_Y: List[float],
        AHE: List[float],
    ) -> None:
        """Creates a junction with a STRIP magnetoresistance.
        Each of the Rx0, Ry, AMR, AMR and SMR is list matching the
        length of the layers passed (they directly correspond to each layer).
        Calculates the magnetoresistance as per: __see reference__:
        Spin Hall magnetoresistance in metallic bilayers by Kim, J. et al.
        :param Rx0
        :param Ry0
        :param AMR_X
        :param AMR_Y
        :param SMR_X
        :param SMR_Y
        :param AHE
        """

    @overload
    def __init__(*args, **kwargs) -> Any:
        ...

    def clearLog(self) -> Dict[str, Any]:
        """
        Reset current simulation state`
        """
        ...

    def getLayerMagnetisation(self, layer_id: str) -> CVector:
        ...

    def getLog(self) -> Dict[str, List[float]]:
        """
        Retrieve the simulation log [data].
        """
        ...

    def getMagnetoresistance(self) -> List[float]:
        ...

    def runSimulation(
        self,
        totalTime: float,
        timeStep: float = ...,
        writeFrequency: float = ...,
        persist: bool = ...,
        log: bool = ...,
        calculateEnergies: bool = ...,
    ) -> None:
        """Main run simulation function.
        Use it to run the simulation.
        :param totalTime: total time of a simulation, give it in seconds. Typical length is in ~couple ns.
        :param timeStep: the integration step of the RK45 method. Default is 1e-13
        :param writeFrequency: how often is the log saved to? Must be no smaller than `timeStep`. Default is 1e-11.
        :param persist: whether to save to the filename specified in the Junction constructor. Default is true
        :param log: if you want some verbosity like timing the simulation. Default is false
        :param calculateEnergies: [WORK IN PROGRESS] log energy values to the log. Default is false.
        """
        ...

    def setIECDriver(self, bottom_layer: str, top_layer: str,
                     driver: ScalarDriver) -> None:
        """Set IEC interaction between two layers.
        The names of the params are only for convention. The IEC will be set
        between bottomLyaer or topLayer, order is irrelevant.
        :param bottomLayer: the first layer id
        :param topLayer: the second layer id
        """
        ...

    def setQuadIECDriver(self, bottom_layer: str, top_layer: str,
                         driver: ScalarDriver) -> None:
        """Set secondary (biquadratic term) IEC interaction between two layers.
        The names of the params are only for convention. The IEC will be set
        between bottomLyaer or topLayer, order is irrelevant.
        :param bottomLayer: the first layer id
        :param topLayer: the second layer id
        """
        ...

    def setLayerTemperatureDriver(self, layer_id: str,
                                  driver: ScalarDriver) -> None:
        ...

    def setLayerAnisotropyDriver(self, layer_id: str,
                                 driver: ScalarDriver) -> None:
        ...

    def setLayerCurrentDriver(self, layer_id: str,
                              driver: ScalarDriver) -> None:
        ...

    def setLayerExternalFieldDriver(self, layer_id: str,
                                    driver: AxialDriver) -> None:
        ...

    def setLayerMagnetisation(self, layer_id: str, mag: CVector) -> None:
        ...

    @overload
    def setLayerOerstedFieldDriver(self, layer_id: str,
                                   driver: AxialDriver) -> None:
        ...

    def setLayerDampingLikeTorqueDriver(self, layer_id: str,
                                        driver: ScalarDriver) -> None:
        """Set the damping like torque driver for a layer.
        :param layer_id: the layer id
        :param driver: the driver
        """
        ...

    def setLayerFieldLikeTorqueDriver(self, layer_id: str,
                                      driver: ScalarDriver) -> None:
        """Set the field like torque driver for a layer.
        :param layer_id: the layer id
        :param driver: the driver
        """
        ...

    def setLayerOneFNoise(self, layer_id: str, sources: int, bias: float,
                          scale: float) -> None:
        """Set 1/f noise for a layer.
        :param layer_id: the layer id
        :param sources: the number of generation sources (the more the slower, but more acc.)
        :param bias: the bias of the noise (p in the Multinomial distribution)
        :param scale: the scale of the noise, additional scaling factor
        """
        ...


class Layer:
    def __init__(
        self,
        id: str,
        mag: CVector,
        anis: CVector,
        Ms: float,
        thickness: float,
        cellSurface: float,
        demagTensor: List[CVector],
        temperature: float = ...,
        damping: float = ...,
    ) -> None:
        """
        The basic structure is a magnetic layer.
        Its parameters are defined by the constructor and may be altered
        by the drivers during the simulation time.
        If you want STT, remember to set the reference vector for the polarisation of the layer.
        Use `setReferenceLayer` function to do that.
        :param id: identifiable name for a layer -- e.g. "bottom" or "free".
        :param mag: initial magnetisation. Must be normalised (norm of 1). Used for quicker convergence.
        :param anis: anisotropy of the layer. A normalised vector
        :param Ms: magnetisation saturation. Unit: Tesla [T].
        :param thickness: thickness of the layer. Unit: meter [m].
        :param cellSurface: surface of the layer, for volume calculation. Unit: meter^2 [m^2].
        :param dipoleTensor: dipole tensor of the layer.
        :param damping: often marked as alpha in the LLG equation. Damping of the layer. Default 0.011. Dimensionless
        """
        ...

    @staticmethod
    def createSOTLayer(
        id: str,
        mag: CVector,
        anis: CVector,
        Ms: float,
        thickness: float,
        cellSurface: float,
        demagTensor: List[CVector],
        damping: float = 0.11,
        fieldLikeTorque: float = 0,
        dampingLikeTorque: float = 0,
    ) -> "Layer":
        """
        Create SOT layer -- including damping and field-like torques that are
        calculated based on the effective Spin Hall angles.
        :param id: identifiable name for a layer -- e.g. "bottom" or "free".
        :param mag: initial magnetisation. Must be normalised (norm of 1). Used for quicker convergence.
        :param anis: anisotropy of the layer. A normalised vector
        :param Ms: magnetisation saturation. Unit: Tesla [T].
        :param thickness: thickness of the layer. Unit: meter [m].
        :param cellSurface: surface of the layer, for volume calculation. Unit: meter^2 [m^2].
        :param dipoleTensor: dipole tensor of the layer.
        :param temperature: resting temperature of the layer. Unit: Kelvin [K].
        :param damping: often marked as alpha in the LLG equation. Damping of the layer. Default 0.011. Dimensionless.
        """
        ...

    @staticmethod
    def createSTTLayer(
        id: str,
        mag: CVector,
        anis: CVector,
        Ms: float,
        thickness: float,
        cellSurface: float,
        demagTensor: List[CVector],
        damping: float = 0.011,
        SlonczewskiSpacerLayerParameter: float = 1.0,
        beta: float = 0.0,
        spinPolarisation: float = 0.0,
    ) -> "Layer":
        """
        Create STT layer -- with the standard Slomczewski formulation.
        :param id: identifiable name for a layer -- e.g. "bottom" or "free".
        :param mag: initial magnetisation. Must be normalised (norm of 1). Used for quicker convergence.
        :param anis: anisotropy of the layer. A normalised vector
        :param Ms: magnetisation saturation. Unit: Tesla [T].
        :param thickness: thickness of the layer. Unit: meter [m].
        :param cellSurface: surface of the layer, for volume calculation. Unit: meter^2 [m^2].
        :param dipoleTensor: dipole tensor of the layer.
        :param damping: often marked as alpha in the LLG equation. Damping of the layer. Default 0.011. Dimensionless.
        :param SlonczewskiSpacerLayerParameter: Slomczewski parameter. Often marked as lambda.
        :param beta: beta parameter that scales FL/DL ratio.
        :param spinPolarisation: the spin effectiveness.
        """
        ...

    def setAnisotropyDriver(self, driver: ScalarDriver) -> None:
        """Set anisotropy driver for the layer.
        It's scalar. The axis is determined in the layer constructor"""
        ...

    def setTemperatureDriver(self, driver: ScalarDriver) -> None:
        """Set a driver for the temperature of the layer.
        Automatically changes the solver to Euler-Heun."""
        ...

    def setExternalFieldDriver(self, driver: AxialDriver) -> None:
        ...

    def setMagnetisation(self, mag: CVector) -> None:
        ...

    def setOerstedFieldDriver(self, driver: AxialDriver) -> None:
        ...

    def setDampingLikeTorqueDriver(self, arg0: ScalarDriver) -> None:
        """Set a driver for the damping like torque of the layer."""
        ...

    def setFieldLikeTorqueDriver(self, arg0: ScalarDriver) -> None:
        """Set a driver for the field like torque of the layer."""
        ...

    def setReferenceLayer(self, ref: CVector) -> None:
        ...

    @overload
    def setReferenceLayer(self, ref: "Reference") -> None:
        ...

    def setTopDipoleTensor(self, tensor: List[CVector]) -> None:
        """Set a dipole tensor from the top layer."""
        ...

    def setBottomDipoleTensor(self, tensor: List[CVector]) -> None:
        """Set a dipole tensor from the bottom layer."""
        ...

    def getId(self) -> str:
        """Get Id of the layer"""
        ...

    def setAlternativeSTT(self, arg0: bool) -> None:
        """Switch to an alternative STT forumulation (Taniguchi et al.)
        https://iopscience.iop.org/article/10.7567/APEX.11.013005
        """
        ...


class NullDriver(ScalarDriver):
    def __init__(self) -> None:
        """
        An empty driver that does nothing. Use in Axial Driver when
        the axis is to be id.
        """
        ...


class ScalarDriver:
    def __init__(self, *args, **kwargs) -> None:
        ...

    @staticmethod
    def getConstantDriver(constantValue: float) -> "ScalarDriver":
        """
        Constant driver produces a constant signal of a fixed amplitude.
        :param constantValue: constant value of the driver (constant offset/amplitude)

        """
        ...

    @staticmethod
    def getPulseDriver(constantValue: float, amplitude: "ScalarDriver",
                       period: float, cycle: float) -> Any:
        """
        Produces a square pulse of certain period and cycle
        :param constantValue: offset (vertical) of the pulse. The pulse amplitude will be added to this.
        :param amplitude: amplitude of the pulse signal
        :param period: period of the signal in seconds
        :param cycle: duty cycle of the signal -- a fraction between [0 and 1].

        """
        ...

    @staticmethod
    def getSineDriver(constantValue: float, amplitude: "ScalarDriver",
                      frequency: float, phase: float) -> Any:
        """
        Produces a sinusoidal signal with some offset (constantValue), amplitude frequency and phase offset.
        :param constantValue: vertical offset. The sine will oscillate around this value.
        :param amplitude: amplitude of the sine wave
        :param frequency: frequency of the sine
        :param phase: phase of the sine in radians.

        """
        ...

    @staticmethod
    def getStepDriver(constantValue: float, amplitude: float, timeStart: float,
                      timeStop: float) -> "ScalarDriver":
        """
        Get a step driver. It has amplitude between timeStart and timeStop and 0 elsewhere
        :param constantValue: offset of the pulse (vertical)
        :param amplitude: amplitude that is added on top of the constantValue
        :param timeStart: start of the pulse
        :param timeStop: when the pulse ends
        """
        ...

    def getTrapezoidDriver(self, constantValue: float, amplitude: float,
                           timeStart, edgeTime: float,
                           steadyTime: float) -> Any:
        """Create Trapezoid driver. Has a rising and a falling edge.
        :param constantValue: offset of the pulse (vertical)
        :param amplitude: amplitude that is added on top of the constantValue
        :param timeStart: start of the pulse
        :param edgeTime: time it takes to reach the maximum amplitude
        :param steadyTime: time it spends in a steady state
        """
        ...


class SolverMode:
    """SolverMode Indicator"""
    DormandPrice: ClassVar[SolverMode] = ...
    EulerHeun: ClassVar[SolverMode] = ...
    RK4: ClassVar[SolverMode] = ...
    Heun: ClassVar[SolverMode] = ...

class Reference:
    """Reference layer indicator."""
    bottom: ClassVar[Reference] = ...
    fixed: ClassVar[Reference] = ...
    none: ClassVar[Reference] = ...
    top: ClassVar[Reference] = ...


DormandPrice: SolverMode
EulerHeun: SolverMode
Heun: SolverMode
RK4: SolverMode
bottom: Reference
fixed: Reference
none: Reference
top: Reference
