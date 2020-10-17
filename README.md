# CMTJ

CMTJ is a set of C++ headers for spin magnetics in Magnetic Tunnel Junctions, with some bindings from C++ to Python, but there's also a Python standalone library called PyMTJ. The Python standalone library, the PyMTJ, is not maintained actively, but you are welcome to contribute if you would like to see it revived. Until then, C++ with Python bindings, i.e. CMTJ, is the best choice if you want to go with a maintained version.


## Installation steps

The easiest way to compile the PyBind11 bindings is to use the provided [Dockerfile](Dockerfile), if you know how to use it:
```bash
# docker build
docker build -t pymtj .
# docker run 
docker run -it --entrypoint bash cmtj 
```
You may want to develop C++ code inside the Dockerfile since it has everything installed (including the `vim`). Then, you may want to mount the volume to container's `/mnt`:
```bash 
docker run -it --entrypoint bash -v /path/to/my/volume:/mnt cmtj
```

### Requirements 
The C++ version basically requires only:  
* FFTW3  
  The library for computing the Fast Fourier Transform, needed for some of the frequency analysis included in the package

  The library may be downloaded from the link here [http://www.fftw.org/](http://www.fftw.org/).  
  However, if you're on the Ubuntu (it may work for some other distros as well), you may install the package with:
  ```bash
  apt-get install -y libfftw3-dev
  ```

* PyBind11 
  PyBind11 allows us to bind pieces of C++ code to be usable from within the Python code.
  This package *is not* necessary if you just plan to use the C++ API. Otherwise, if you  need to compile the PyMTJ yourself, then you should install this package.

  The library is available here [https://github.com/pybind/pybind11](https://github.com/pybind/pybind11).


### Compilation
_The header files are available in the [`cmtj`](cmtj/) folder_

The actual installation is rather simple -- given the Makefile I have provided, you may just compile your code using one of the MakeFile commands since the C++ API is virtually entirely header (2 headers -- [`junction.hpp`](cmtj/junction.hpp) and [`cvector.hpp`](cmtj/cvector.hpp). Simply include the headers in your code and compile it away.

The only actual compilation is required for the Python bindings. You may want to use one of the `make python` command variations from the Makefile.

_Which one to run it?_   
For the MacOs use `make python-macos`, for Ubuntu use `make python-ubuntu`

### Easy binding installation 
The bindings to Python are be installed into the Python by running the:
```bash 
pip3 install -e .
```
in the [`cmtj`](cmtj/) folder of the repository -- it will automatically install the cmtj module into Python.   
Alternatively, you may run 
```bash 
python3 setup.py build_ext --inplace
```
in the same, [`cmtj`](cmtj/), folder.

Now, you may use the Python bindings as any other module in Python. See the [examples] (`examples/`) folder for some instruction on how to use the module.


------------------
# PYTMJ 

## Python version (Python libray may be found in [`pymtj`](pymtj/) folder)
Below is a short introduction to that library, which is purely Python (no bindings, simply Cython).
A Python library for spin magnetics in Magnetic Tunnel Junctions.

### Installation
Installation is simple and requires Python3 (preferably 3.7)
Clone this repository
```bash 
git clone https://github.com/LemurPwned/spinpy.git
```
Enter the [`pymtj`](pymtj/) repository and: 
```bash 
pip3 install -e .
```
This will install and the dependencies automatically and all the dependencies as well.
For the precise requirements that will be installed please see `setup.py` files and the 
requirements section.


### Using the software
Please look into the [`pymtj/automation/workers.py`](pymtj/automation/workers.py) for the advanced examples of usage 
and [`test_run.py`](test_run.py) for an example on how to run test automated tasks

Please fill the issues or pull requests in case you require any additional functionality 
or encouter any problems with the software
#### Warning

Several functions allow you to set an external
Python function that updates the state of the junction
such as `.set_global_coupling_function` or `set_junction_global_external_field`  
Be very careful because setting the improper values might cause the prolonged computation times!
This library is optimized with Cython to serve the Python's `float` object which is 24 bytes. What it means for computation is that when you pass it a 32 byte object such as `numpy.float64` or even `numpy.float32`, the calculation will take significantly much more!
Take a note that for accurate results in macromagnetics you do not need such a large precision, thus you might be well off with this precision loss, expecially for the external field values.

### Coupling 

Setting the coupling value in the Junction `couplings` is defaulted to `[[2], [1]]`. Note, that the coupling is always respective to the layer ids. Layer ids start from 1 to N, where N is the last layer. Hence, `[[2], [1]]` means that the first layer is coupled with `2` layer and second layer is coupled with `1` layer.
Otherwise it can be read:
```python
[    
    [2, 3, 4] # layer 1 couples with 2, 3, 4
    [1, 3] # layer 2 couples with 1, 3,
    ... # etc. (must be provided for all layrs)
]
```

## **Calculations**

### Conventional markings   
1. $\vec{m}_i$      the magnetisation vector of the $i$th layer. If not given, means the *current* layer 
2. $M_s$            is the magnetisation saturation 
3. $\vec{H}$        is the vector field  
### Demagnetisation tensor 
$$\vec{H_i} = -NM_s^{j}\vec{m_i}$$
where the $M_s^{j}$ means the saturation of the *other*, *non-local* layer $j$ and $N$ is the demagnetisation tensor dependent on the elative shape of the layers $i$ and $j$.

### Interlayer Exchange Coupling interaction (IEC)

$$\vec{H_i} = \frac{J}{\mu M_s d} (\vec{m_j} - \vec{m_i})$$
where the $J$ is the coupling constant between the layers $i$ and $j$, and $d$ is the thickness.
