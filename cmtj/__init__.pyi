import typing
from typing import Any, ClassVar, overload

xaxis: Axis
yaxis: Axis
zaxis: Axis
none: Axis
all: Axis

def c_dot(arg0: CVector, arg1: CVector) -> float:
    """Compute dot (scalar) product of two CVectors."""
    ...

def constantDriver(constant: float) -> ScalarDriver:
    """
    Constant driver produces a constant signal of a fixed amplitude.
    :param constant: constant value of the driver (constant offset/amplitude)
    """
    ...

def sineDriver(constantValue: float, amplitude: float, frequency: float, phase: float) -> ScalarDriver:
    """
    Produces a sinusoidal signal with some offset (constantValue), amplitude frequency and phase offset.
    :param constantValue: vertical offset. The sine will oscillate around this value.
    :param amplitude: amplitude of the sine wave
    :param frequency: frequency of the sine
    :param phase: phase of the sine in radians.
    """
    ...

def gaussianImpulseDriver(constantValue: float, amplitude: float, t0: float, sigma: float) -> ScalarDriver:
    """
    Gaussian impulse driver. It starts with an max amplitude at t0 and falls off with sigma.

    Formula:

    $A \exp(-(t - t_0)^2 / (2\sigma^2))$

    :param constantValue: offset of the pulse (vertical)
    :param amplitude: amplitude that is added on top of the constantValue
    :param t0: start of the pulse
    :param sigma: fall-off of the Gaussian pulse
    """
    ...

def gaussianStepDriver(constantValue: float, amplitude: float, t0: float, sigma: float) -> ScalarDriver:
    """Gaussian step driver (erf function). It starts at t0 and falls off with sigma.

    Formula:

    $f(t) = c + A + A\mathrm{erf}((t - t_0) / (\sigma \sqrt(2)))$

    :param constantValue: offset of the pulse (vertical)
    :param amplitude: amplitude that is added on top of the constantValue
    :param t0: start of the pulse
    :param sigma: fall-off of the Gaussian pulse
    """
    ...

def posSineDriver(constantValue: float, amplitude: float, frequency: float, phase: float) -> ScalarDriver:
    """Produces a positive sinusoidal signal with some offset (constantValue), amplitude frequency and phase offset.
    :param constantValue: vertical offset. The sine will oscillate around this value.
    :param amplitude: amplitude of the sine wave
    :param frequency: frequency of the sine
    :param phase: phase of the sine in radians.
    """
    ...

def pulseDriver(constantValue: float, amplitude: float, period: float, cycle: float) -> ScalarDriver:
    """
    Produces a square pulse of certain period and cycle
    :param constantValue: offset (vertical) of the pulse. The pulse amplitude will be added to this.
    :param amplitude: amplitude of the pulse signal
    :param period: period of the signal in seconds
    :param cycle: duty cycle of the signal -- a fraction between [0 and 1].
    """
    ...

def stepDriver(constantValue: float, amplitude: float, timeStart: float, timeStop: float) -> ScalarDriver:
    """
    Get a step driver. It has amplitude between timeStart and timeStop and 0 elsewhere
    :param constantValue: offset of the pulse (vertical)
    :param amplitude: amplitude that is added on top of the constantValue
    :param timeStart: start of the pulse
    :param timeStop: when the pulse ends
    """
    ...

def trapezoidDriver(
    constantValue: float,
    amplitude: float,
    timeStart,
    edgeTime: float,
    steadyTime: float,
) -> ScalarDriver:
    """Create Trapezoid driver. Has a rising and a falling edge.
    :param constantValue: offset of the pulse (vertical)
    :param amplitude: amplitude that is added on top of the constantValue
    :param timeStart: start of the pulse
    :param edgeTime: time it takes to reach the maximum amplitude
    :param steadyTime: time it spends in a steady state
    """
    ...

class AxialDriver:
    """Axial driver class."""

    @overload
    def __init__(self, x: ScalarDriver, y: ScalarDriver, z: ScalarDriver) -> None:
        """Create an axial driver with three scalar drivers for each axis.
        :param x: driver for the x axis
        :param y: driver for the y axis
        :param z: driver for the z axis
        """
        ...

    @overload
    def __init__(self, axialDrivers: list[ScalarDriver]) -> None:
        """Create an axial driver with a list of scalar drivers.
        :param axialDrivers: list of scalar drivers
        """
        ...

    @overload
    def __init__(self, x: float, y: float, z: float) -> None:
        """Create an axial driver with a list of floats.
        :param x: constant float for the x axis
        :param y: constant float for the y axis
        :param z: constant float for the z axis
        """
        ...

    @overload
    def __init__(self, xyz: CVector) -> None:
        """Create an axial driver with a vector.
        :param xyz: CVector object with x, y, z components
        """
        ...

    @overload
    def __init__(*args, **kwargs) -> Any: ...
    @overload
    def applyMask(self, mask: CVector) -> None:
        """Apply mask to the driver.
        :param mask: mask to be applied"""
        ...

    @overload
    def applyMask(self, mask: list[int]) -> None:
        """Apply mask to the driver.
        :param mask: mask to be applied"""
        ...

    @overload
    def applyMask(*args, **kwargs) -> Any: ...
    def getCurrentAxialDrivers(self, arg0: float) -> CVector: ...
    def getVectorAxialDriver(self, arg0: float, arg1: float) -> AxialDriver: ...

class Axis:
    __entries: Any = ...
    xaxis: Any = ...
    yaxis: Any = ...
    zaxis: Any = ...
    none: Any = ...
    all: Any = ...

    def __init__(self, value: int) -> None: ...
    def __eq__(self, other: object) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: object) -> bool: ...
    def __setstate__(self, state: int) -> None: ...
    @property
    def name(self) -> Any: ...
    @property
    def __doc__(self) -> Any: ...
    @property
    def __members__(self) -> Any: ...

class CVector:
    """CVector class. Represents a 3D vector."""

    def __init__(self, x: float, y: float, z: float) -> None:
        """Initialises a 3D vector.
        :param x: x component of the vector
        :param y: y component of the vector
        :param z: z component of the vector"""
        ...

    def length(self) -> float:
        """Returns the length of the vector."""
        ...

    def normalize(self) -> None:
        """Normalizes the vector."""
        ...

    def tolist(self) -> list[float]:
        """Converts the vector to a list."""
        ...

    def __add__(self, arg0: CVector) -> CVector: ...
    def __eq__(self, arg0: CVector) -> bool: ...
    def __getitem__(self, arg0: int) -> float: ...
    def __iter__(self) -> typing.Iterator[float]: ...
    def __iadd__(self, arg0: CVector) -> CVector: ...
    def __imul__(self, arg0: float) -> CVector: ...
    def __isub__(self, arg0: CVector) -> CVector: ...
    def __len__(self) -> int: ...
    def __mul__(self, arg0: float) -> CVector: ...
    def __ne__(self, arg0: CVector) -> bool: ...
    def __rmul__(self, arg0: float) -> CVector: ...
    def __sub__(self, arg0: CVector) -> CVector: ...
    @property
    def x(self) -> float: ...
    @x.setter
    def x(self, val: float) -> None: ...
    @property
    def y(self) -> float: ...
    @y.setter
    def y(self, val: float) -> None: ...
    @property
    def z(self) -> float: ...
    @z.setter
    def z(self, val: float) -> None: ...

class Junction:
    @overload
    def __init__(self, layers: list[Layer]) -> None:
        """"""
        ...

    @overload
    def __init__(self, layers: list[Layer], Rp: float = ..., Rap: float = ...) -> None:
        """Creates a junction with a magnetoresistance.
        :param layers: list of layers

        :param Rp: Parallel magnetoresistance
        :param Rap: Magnetoresistance anti-parallel state
        """
        ...

    @overload
    def __init__(
        self,
        layers: list[Layer],
        Rx0: list[float],
        Ry0: list[float],
        AMR_X: list[float],
        AMR_Y: list[float],
        SMR_X: list[float],
        SMR_Y: list[float],
        AHE: list[float],
    ) -> None:
        """Creates a junction with a STRIP magnetoresistance.
        Each of the Rx0, Ry, AMR, AMR and SMR is list matching the
        length of the layers passed (they directly correspond to each layer).
        Calculates the magnetoresistance as per: __see reference__:
        Spin Hall magnetoresistance in metallic bilayers by Kim, J. et al.
        :param Rx0: Magnetoresistance offset longitudinal
        :param Ry0: Magnetoresistance offset transverse
        :param AMR_X: Anisotropic magnetoresistance longitudinal
        :param AMR_Y: Anisotropic magnetoresistance transverse
        :param SMR_X: Spin magnetoresistance longitudinal
        :param SMR_Y: Spin magnetoresistance transverse
        :param AHE: Anomalous Hall effect resistance offset (transverse only)
        """
        ...

    @overload
    def __init__(*args, **kwargs) -> Any: ...
    def clearLog(self) -> dict[str, Any]:
        """Reset current simulation state."""
        ...

    def getLayerMagnetisation(self, layerId: str) -> CVector:
        """Get the magnetisation of a layer.
        :param layerId: the layer id"""
        ...

    def getLog(self) -> dict[str, list[float]]:
        """Retrieve the simulation log [data]."""
        ...

    def getMagnetoresistance(self) -> list[float]: ...
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

    def setIECDriver(self, bottomLayer: str, topLayer: str, driver: ScalarDriver) -> None:
        """Set IEC interaction between two layers.
        The names of the params are only for convention. The IEC will be set
        between bottomLyaer or topLayer, order is irrelevant.
        :param bottomLayer: the first layer id
        :param topLayer: the second layer id
        """
        ...

    def setQuadIECDriver(self, bottomLayer: str, topLayer: str, driver: ScalarDriver) -> None:
        """Set secondary (biquadratic term) IEC interaction between two layers.
        The names of the params are only for convention. The IEC will be set
        between bottomLyaer or topLayer, order is irrelevant.
        :param bottomLayer: the first layer id
        :param topLayer: the second layer id
        """
        ...

    def setLayerTemperatureDriver(self, layerId: str, driver: ScalarDriver) -> None:
        """Set a temperature driver for a layer.
        :param layerId: the id of the layer.
        :param driver: the temperature driver to be set.
        """
        ...

    def setLayerAnisotropyDriver(self, layerId: str, driver: ScalarDriver) -> None:
        """Set anisotropy driver for a layer.
        :param layerId: the id of the layer.
        :param driver: the anisotropy driver to be set.
        """
        ...

    def setLayerCurrentDriver(self, layerId: str, driver: ScalarDriver) -> None:
        """Set a current driver for a layer.
        :param layerId: the layer id
        :param driver: the driver
        """
        ...

    def setLayerExternalFieldDriver(self, layerId: str, driver: AxialDriver) -> None:
        """Set an external field driver for a layer.
        :param layerId: the id of the layer.
        :param driver: the field driver to be set.
        """
        ...

    def setLayerMagnetisation(self, layerId: str, mag: CVector) -> None:
        """Set the magnetisation of a layer.
        :param layerId: the layer id
        :param mag: the magnetisation
        """
        ...

    @overload
    def setLayerOerstedFieldDriver(self, layerId: str, driver: AxialDriver) -> None:
        """Set an Oersted field driver for a layer.
        :param layerId: the id of the layer.
        :param driver: the field driver to be set.
        """
        ...

    def setLayerDampingLikeTorqueDriver(self, layerId: str, driver: ScalarDriver) -> None:
        """Set the damping like torque driver for a layer.
        :param layerId: the layer id
        :param driver: the driver
        """
        ...

    def setLayerFieldLikeTorqueDriver(self, layerId: str, driver: ScalarDriver) -> None:
        """Set the field like torque driver for a layer.
        :param layerId: the layer id
        :param driver: the driver
        """
        ...

    def setLayerOneFNoise(self, layerId: str, sources: int, bias: float, scale: float) -> None:
        """Set 1/f noise for a layer.
        :param layerId: the layer id
        :param sources: the number of generation sources (the more the slower, but more acc.)
        :param bias: the bias of the noise (p in the Multinomial distribution)
        :param scale: the scale of the noise, additional scaling factor
        """
        ...

    def getLayer(self, layerId: str) -> Layer:
        """Get a specific layer from the junction. Returns a reference.
        :param layerId: the id of the layer (string) as passed in the init.
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
        demagTensor: list[CVector],
        temperature: float = ...,
        damping: float = ...,
    ) -> Layer:
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
        demagTensor: list[CVector],
        damping: float = 0.11,
        fieldLikeTorque: float = 0,
        dampingLikeTorque: float = 0,
    ) -> Layer:
        """
        Create SOT layer -- including damping and field-like torques that are
        calculated based on the effective Spin Hall angles.
        :param id: identifiable name for a layer -- e.g. "bottom" or "free".
        :param mag: initial magnetisation. Must be normalised (norm of 1). Used for quicker convergence.
        :param anis: anisotropy of the layer. A normalised vector
        :param Ms: magnetisation saturation. Unit: Tesla [T].
        :param thickness: thickness of the layer. Unit: meter [m].
        :param cellSurface: surface of the layer, for volume calculation. Unit: meter^2 [m^2].
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
        demagTensor: list[CVector],
        damping: float = 0.011,
        SlonczewskiSpacerLayerParameter: float = 1.0,
        beta: float = 0.0,
        spinPolarisation: float = 0.0,
    ) -> Layer:
        """
        Create STT layer -- with the standard Slomczewski formulation.
        :param id: identifiable name for a layer -- e.g. "bottom" or "free".
        :param mag: initial magnetisation. Must be normalised (norm of 1). Used for quicker convergence.
        :param anis: anisotropy of the layer. A normalised vector
        :param Ms: magnetisation saturation. Unit: Tesla [T].
        :param thickness: thickness of the layer. Unit: meter [m].
        :param cellSurface: surface of the layer, for volume calculation. Unit: meter^2 [m^2].
        :param damping: often marked as alpha in the LLG equation. Damping of the layer. Default 0.011. Dimensionless.
        :param SlonczewskiSpacerLayerParameter: Slomczewski parameter. Often marked as lambda.
        :param beta: beta parameter that scales FL/DL ratio.
        :param spinPolarisation: the spin effectiveness.
        """
        ...

    def createBufferedAlphaNoise(self, bufferSize: int) -> None:
        """Create a buffered alpha noise generator."""
        ...

    def setAlphaNoise(self, alpha: float, std: float, scale: float, axis: Axis = Axis.all) -> None:
        """Set alpha noise for the layer.
        :param alpha: Alpha parameter
        :param std: Standard deviation
        :param scale: Scale
        :param axis: Axis, by default all axes are used
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

    def setExternalFieldDriver(self, driver: AxialDriver) -> None: ...
    def setMagnetisation(self, mag: CVector) -> None:
        """Set the magnetisation of the layer.
        :param mag: the magnetisation to be set."""
        ...

    def setOerstedFieldDriver(self, driver: AxialDriver) -> None:
        """Set an Oersted field driver for the layer.
        :param driver: the field driver to be set."""
        ...

    def setDampingLikeTorqueDriver(self, driver: ScalarDriver) -> None:
        """Set a driver for the damping like torque of the layer.
        :param driver: the driver to be set."""
        ...

    def setFieldLikeTorqueDriver(self, driver: ScalarDriver) -> None:
        """Set a driver for the field like torque of the layer.
        :param driver: the driver to be set."""
        ...

    def setReferenceLayer(self, ref: CVector) -> None:
        """Set a reference layer for the STT.
        :param ref: the reference layer vector."""
        ...

    @overload
    def setReferenceLayer(self, ref: Reference) -> None:  # noqa: F811
        """Set a reference layer for the STT. The reference can be
        FIXED, BOTTOM or TOP. YOu can use another layer as reference
        to this one.
        :param ref: the reference layer vector."""
        ...

    def setTopDipoleTensor(self, tensor: list[CVector]) -> None:
        """Set a dipole tensor from the top layer.
        :param tensor: the dipole tensor to be set.
        """
        ...

    def setBottomDipoleTensor(self, tensor: list[CVector]) -> None:
        """Set a dipole tensor from the bottom layer.
        :param tensor: the dipole tensor to be set.
        """
        ...

    def getId(self) -> str:
        """Get Id of the layer"""
        ...

    def setAlternativeSTT(self, setAlternative: bool) -> None:
        """Switch to an alternative STT forumulation (Taniguchi et al.)
        https://iopscience.iop.org/article/10.7567/APEX.11.013005
        :param setAlternative: whether to set the alternative STT formulation
        """
        ...

    def setKappa(self, kappa: float) -> None:
        """Set the kappa parameter for the layer -- determines SOT mixing
            Hdl * kappa + Hfl
        Allows you to turn off Hdl. Turning Hfl is via beta parameter.
        :param kappa: the kappa parameter
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
    def __init__(self, *args, **kwargs) -> None: ...
    def getCurrentScalarValue(self, time: float) -> float:
        """
        :param time: time in seconds
        :return: the scalar value of the driver at time.
        """
        ...

    @staticmethod
    def getConstantDriver(constantValue: float) -> ScalarDriver:
        """
        Constant driver produces a constant signal of a fixed amplitude.
        :param constantValue: constant value of the driver (constant offset/amplitude)
        """
        ...

    @staticmethod
    def getPulseDriver(constantValue: float, amplitude: float, period: float, cycle: float) -> ScalarDriver:
        """
        Produces a square pulse of certain period and cycle
        :param constantValue: offset (vertical) of the pulse. The pulse amplitude will be added to this.
        :param amplitude: amplitude of the pulse signal
        :param period: period of the signal in seconds
        :param cycle: duty cycle of the signal -- a fraction between [0 and 1].
        """
        ...

    @staticmethod
    def getSineDriver(constantValue: float, amplitude: ScalarDriver, frequency: float, phase: float) -> Any:
        """
        Produces a sinusoidal signal with some offset (constantValue), amplitude frequency and phase offset.
        :param constantValue: vertical offset. The sine will oscillate around this value.
        :param amplitude: amplitude of the sine wave
        :param frequency: frequency of the sine
        :param phase: phase of the sine in radians.
        """
        ...

    @staticmethod
    def getStepDriver(constantValue: float, amplitude: float, timeStart: float, timeStop: float) -> ScalarDriver:
        """
        Get a step driver. It has amplitude between timeStart and timeStop and 0 elsewhere
        :param constantValue: offset of the pulse (vertical)
        :param amplitude: amplitude that is added on top of the constantValue
        :param timeStart: start of the pulse
        :param timeStop: when the pulse ends
        """
        ...

    @staticmethod
    def getTrapezoidDriver(
        constantValue: float,
        amplitude: float,
        timeStart,
        edgeTime: float,
        steadyTime: float,
    ) -> ScalarDriver:
        """Create Trapezoid driver. Has a rising and a falling edge.
        :param constantValue: offset of the pulse (vertical)
        :param amplitude: amplitude that is added on top of the constantValue
        :param timeStart: start of the pulse
        :param edgeTime: time it takes to reach the maximum amplitude
        :param steadyTime: time it spends in a steady state
        """
        ...

    @staticmethod
    def getGaussianImpulseDriver(constantValue: float, amplitude: float, t0: float, sigma: float) -> ScalarDriver:
        """Gaussian impulse driver. It has amplitude starts at t0 and falls off with sigma.

        Formula:
        A * exp(-((t - t0) ** 2) / (2 * sigma ** 2))

        :param constantValue: offset of the pulse (vertical)
        :param amplitude: amplitude that is added on top of the constantValue
        :param t0: start of the pulse
        :param sigma: fall-off of the Gaussian pulse
        """
        ...

    @staticmethod
    def getGaussianStepDriver(constantValue: float, amplitude: float, t0: float, sigma: float) -> ScalarDriver:
        """Gaussian step driver (erf function). It has amplitude starts at t0 and falls off with sigma.

        Formula:
        f(t) = constantValue + amplitude * (1 + erf((t - t0) / (sigma * sqrt(2))))

        :param constantValue: offset of the pulse (vertical)
        :param amplitude: amplitude that is added on top of the constantValue
        :param t0: start of the pulse
        :param sigma: fall-off of the Gaussian pulse
        """
        ...

    @staticmethod
    def getPosSineDriver(constantValue: float, amplitude: float, frequency: float, phase: float) -> ScalarDriver:
        """Produces a positive sinusoidal signal with some offset (constantValue), amplitude frequency and phase offset.
        :param constantValue: vertical offset. The sine will oscillate around this value.
        :param amplitude: amplitude of the sine wave
        :param frequency: frequency of the sine
        :param phase: phase of the sine in radians.
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
