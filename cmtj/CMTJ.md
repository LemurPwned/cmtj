# CMTJ
## Brief description
CMTJ is a library for simulating multilayer magnetic structures.
Each structure is called a `junction` and is composed of magnetic `layers`.

## Installation 
Both C++ and Python interfaces are first-class citizens, but because of how this project is structured 
Python code is generated from the C++ code. In case of linting problems, please report the issue on Github.
### C++
The C++ library requires no installation, only the inclusion of the headers.
The basic headers you'll need are [`junction.hpp`](core/junction.hpp) and [`cvector.hpp`](core/cvector.hpp).
Beyond that, there's [`compute.hpp`](core/compute.hpp) that contains utilities for FFT compute, mainly for Voltage Spin Diode and PIMM.
There's also [`stack.hpp`](core/stack.hpp), for stacking the whole junctions into stacks, but as of now it's WIP. 

But if you wish to use functions from [`compute.hpp`](core/compute.hpp), then there's an extra dependency on `fftw3`. 
To install it follow the main source page for `fftw3` -- [http://www.fftw.org/](http://www.fftw.org/).
Here are some quick things that may work for you.
* MacOs
  ```bash
  brew install fftw
  ```
* Linux
  ```bash
  apt-get install -y libfftw3-dev
  ```
* Windows 
  (no idea, please check the fftw source page linked above)

### Python
The python package is to be installed from main directory. Go up from this file's directory. 
Then, launch the pip installation. It should compile the C++ code and bind it to the Python module.
Please report all compilation errors so that we may help you debug it.

You'll also need `pybind11` python package:
```bash 
pip3 install pybind11
``` 

Pip installation command:
```bash
pip3 install -e .
```


## Requirements summary :warning:
The C++ version basically requires only:  
* FFTW3  (only for [`compute.hpp`](core/compute.hpp))
  The library for computing the Fast Fourier Transform, needed for some of the frequency analysis included in the package

  The library may be downloaded from the link here [http://www.fftw.org/](http://www.fftw.org/).  
  However, if you're on the Ubuntu (it may work for some other distros as well), you may install the package with:
  ```bash
  apt-get install -y libfftw3-dev
  ```

* PyBind11 (only if you want to compile the `python3` yourself)
  PyBind11 allows us to bind pieces of C++ code to be usable from within the Python code.
  This package *is not* necessary if you just plan to use the C++ API. Otherwise, if you  need to compile the PyMTJ yourself, then you should install this package.

  The library is available here [https://github.com/pybind/pybind11](https://github.com/pybind/pybind11).


## CMTJ-SERVER Now Available :boom:
For some real-world use cases (and a nice server wrapper too:exclamation:), please check the *[cmtj-server project](https://github.com/LemurPwned/cmtj-server)*.   
It contains some ready-to-use examples for the experimental setting. If you'd like me to add an experimental setup, post an issue there. For theorethical contributions, submit an issue in this repo.
