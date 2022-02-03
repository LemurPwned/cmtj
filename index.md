
<img style="float: right; max-width: 50px;" src="./assets/icon.svg">

[![PyPI](https://github.com/LemurPwned/cmtj/actions/workflows/main.yml/badge.svg?branch=master)](https://github.com/LemurPwned/cmtj/actions/workflows/main.yml)
[![pages-build-deployment](https://github.com/LemurPwned/cmtj/actions/workflows/pages/pages-build-deployment/badge.svg?branch=gh-pages)](https://github.com/LemurPwned/cmtj/actions/workflows/pages/pages-build-deployment)

## Short description
The library allows for macromagnetic simulation of multilayer spintronic structures.
A name may be misleading -- the MTJ (Magnetic Tunnel Junctions) are not the only structures that may be simulated. 

## Quickstart
#### Installation
Installation is as easy as doing:
```bash
python3 -m pip install cmtj
```

If you prefer to clone the repo first and then install directly from the source:
```bash
git clone https://github.com/LemurPwned/cmtj.git
python3 -m pip install .
```
### Examples
Please view the [Examples](#examples) section to get a better grip of the library. Most of the code in there is for the Python bindings but if you prefer, you can always import headers from C++. The function names, library operation will be the same, only the `utils` submodule will be unavailable for the C++ code.  

## Extensions 
There's a GUI version available! If you wish to conduct a subset of simulations, mainly for experimental modelling, please see the *PyMag* project. It uses CMTJ as a backend for fast computation.

## Citing 
Please cite if you decide to use the project
```bibtex 
@article{zietek_numerical_2022,
	title = {Numerical Model Of Harmonic Hall Voltage Detection For Spintronic Devices},
	url = {https://arxiv.org/abs/2202.00364v1},
	author = {Ziętek, Sławomir and Mojsiejuk, Jakub and Grochot, Krzysztof and Łazarski, Stanisław and Skowroński, Witold and Stobiecki, Tomasz},
	urldate = {2022-02-03},
	date = {2022-02-01}
}
```

## Contributions
All contributions are welcome, please leave an issue if you've encountered any trouble with setup or running the library.

