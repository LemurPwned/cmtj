import cmtj

class BufferedAlphaNoise:
    """Create a buffer of alpha noise generator. Alpha can be in [0, 2]."""

    def __init__(self, bufferSize: int, alpha: float, std: float, scale: float) -> None: ...
    def fillBuffer(self) -> None:
        """Fill the buffer with the noise. This method is called only once."""
        ...

    def tick(self) -> float:
        """Produce the next sample of the noise."""
        ...

class VectorAlphaNoise:
    """Create a vector alpha noise generator. Alpha can be in [0, 2]."""

    def __init__(
        self,
        bufferSize: int,
        alpha: float,
        std: float,
        scale: float,
        axis: cmtj.Axis = cmtj.Axis.all,
    ) -> None:
        """Kasdin algorithm for vector alpha noise generation.

        :param bufferSize: Buffer size
        :param alpha: Alpha parameter
        :param std: Standard deviation
        :param scale: Scale
        :param axis: Axis, by default all axes are used
        """
        ...

    def getPrevSample(self) -> cmtj.CVector:
        """Get the previous sample of the noise in a vector form."""
        ...

    def getScale(self) -> float:
        """Get the scale of the noise."""
        ...

    def tick(self) -> float: ...
    def tickVector(self) -> cmtj.CVector:
        """Get the next sample of the noise in a vector form."""
        ...
