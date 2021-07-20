# main module 

```python
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
    def __init__(self, layers: List[Layer], filename: str, Rx0: List[float],
                 Ry0: List[float], AMR_X: List[float], AMR_Y: List[float],
                 SMR_X: List[float], SMR_Y: List[float],
                 AHE: List[float]) -> None:
        """ 
        Creates a junction with a STRIP magnetoresistance.  
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

    def runSimulation(self,
                      totalTime: float,
                      timeStep: float = ...,
                      writeFrequency: float = ...,
                      persist: bool = ...,
                      log: bool = ...,
                      calculateEnergies: bool = ...) -> None:
        """
        Main run simulation function. Use it to run the simulation.
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
        """
        Set IEC interaction between two layers.
        The names of the params are only for convention. The IEC will be set 
        between bottomLyaer or topLayer, order is irrelevant.
        :param bottomLayer: the first layer id
        :param topLayer: the second layer id
        """
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


class Layer:
    def __init__(self,
                 id: str,
                 mag: CVector,
                 anis: CVector,
                 Ms: float,
                 thickness: float,
                 cellSurface: float,
                 demagTensor: List[CVector],
                 dipoleTensor: List[CVector],
                 temperature: float = ...,
                 damping: float = ...) -> None:
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
        :param demagTensor: demagnetisation tensor of the layer.
        :param dipoleTensor: dipole tensor of the layer.
        :param temperature: resting temperature of the layer. Unit: Kelvin [K].
        :param damping: often marked as alpha in the LLG equation. Damping of the layer. Default 0.011. Dimensionless
        """
        ...

    @staticmethod
    def createSOTLayer(id: str,
                       mag: CVector,
                       anis: CVector,
                       Ms: float,
                       thickness: float,
                       cellSurface: float,
                       demagTensor: List[CVector],
                       dipoleTensor: List[CVector],
                       damping: float = 0.11,
                       fieldLikeTorque: float = 0,
                       dampingLikeTorque: float = 0,
                       temperature: float = 0) -> 'Layer':
        """
        Create SOT layer -- including damping and field-like torques that are 
        calculated based on the effective Spin Hall angles.
        :param id: identifiable name for a layer -- e.g. "bottom" or "free".
        :param mag: initial magnetisation. Must be normalised (norm of 1). Used for quicker convergence.
        :param anis: anisotropy of the layer. A normalised vector
        :param Ms: magnetisation saturation. Unit: Tesla [T].
        :param thickness: thickness of the layer. Unit: meter [m].
        :param cellSurface: surface of the layer, for volume calculation. Unit: meter^2 [m^2].
        :param demagTensor: demagnetisation tensor of the layer.
        :param dipoleTensor: dipole tensor of the layer.
        :param temperature: resting temperature of the layer. Unit: Kelvin [K].
        :param damping: often marked as alpha in the LLG equation. Damping of the layer. Default 0.011. Dimensionless.
        """
        ...

    @staticmethod
    def createSTTLayer(id: str,
                       mag: CVector,
                       anis: CVector,
                       Ms: float,
                       thickness: float,
                       cellSurface: float,
                       demagTensor: List[CVector],
                       dipoleTensor: List[CVector],
                       damping: float = 0.011,
                       SlonczewskiSpacerLayerParameter: float = 1.0,
                       beta: float = 0.0,
                       spinPolarisation: float = 0.0,
                       temperature: float = 0) -> 'Layer':
        """
        Create STT layer -- with the standard Slomczewski formulation.
        :param id: identifiable name for a layer -- e.g. "bottom" or "free".
        :param mag: initial magnetisation. Must be normalised (norm of 1). Used for quicker convergence.
        :param anis: anisotropy of the layer. A normalised vector
        :param Ms: magnetisation saturation. Unit: Tesla [T].
        :param thickness: thickness of the layer. Unit: meter [m].
        :param cellSurface: surface of the layer, for volume calculation. Unit: meter^2 [m^2].
        :param demagTensor: demagnetisation tensor of the layer.
        :param dipoleTensor: dipole tensor of the layer.
        :param temperature: resting temperature of the layer. Unit: Kelvin [K].
        :param damping: often marked as alpha in the LLG equation. Damping of the layer. Default 0.011. Dimensionless.
        :param SlonczewskiSpacerLayerParameter: Slomczewski parameter. Often marked as lambda.
        :param beta: beta parameter that scales FL/DL ratio.
        :param spinPolarisation: the spin effectiveness.
        """
        ...

    def setAnisotropyDriver(self, driver: ScalarDriver) -> None:
        ...

    def setExternalFieldDriver(self, driver: AxialDriver) -> None:
        ...

    def setMagnetisation(self, mag: CVector) -> None:
        ...

    def setOerstedFieldDriver(self, driver: AxialDriver) -> None:
        ...

    def setReferenceLayer(self, ref: CVector) -> None:
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
    def getConstantDriver(constantValue: float) -> 'ScalarDriver':
        """
        Constant driver produces a constant signal of a fixed amplitude.
        :param constantValue: constant value of the driver (constant offset/amplitude)

        """
        ...

    @staticmethod
    def getPulseDriver(constantValue: float, amplitude: 'ScalarDriver',
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
    def getSineDriver(constantValue: float, amplitude: 'ScalarDriver',
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
                      timeStop: float) -> 'ScalarDriver':
        """
        Get a step driver. It has amplitude between timeStart and timeStop and 0 elsewhere
        :param constantValue: offset of the pulse (vertical)
        :param amplitude: amplitude that is added on top of the constantValue
        :param timeStart: start of the pulse
        :param timeStop: when the pulse ends
        """
        ...
```