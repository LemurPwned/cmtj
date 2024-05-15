## `BufferedAlphaNoise`

### `__init__(self, bufferSize: int, alpha: float, std: float, scale: float)`

### `fillBuffer(self)`

Fill the buffer with the noise. This method is called only once.

### `tick(self)`

Produce the next sample of the noise.

## `VectorAlphaNoise`

### `__init__(self,bufferSize: int,alpha: float,std: float,scale: float,axis: cmtj.Axis = cmtj.Axis.all,)`

### `getPrevSample(self) -> cmtj.CVector:"""Get the previous sample of the noise in a vector form."""...def getScale(self)`

Get the scale of the noise.

### `tick(self)`
