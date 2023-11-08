<img style="float: right; max-width: 50px;" src="docs/assets/icon.svg">

# CMTJ

[![PyPI](https://github.com/LemurPwned/cmtj/actions/workflows/main.yml/badge.svg?branch=master)](https://github.com/LemurPwned/cmtj/actions/workflows/main.yml)
[![pages-build-deployment](https://github.com/LemurPwned/cmtj/actions/workflows/pages/pages-build-deployment/badge.svg?branch=gh-pages)](https://github.com/LemurPwned/cmtj/actions/workflows/pages/pages-build-deployment)
[![Version](https://img.shields.io/pypi/v/cmtj)](https://pypi.org/project/cmtj/)
[![License](https://img.shields.io/pypi/l/cmtj.svg)](https://github.com/LemurPwned/cmtj/blob/master/LICENSE)
[![Streamlit](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](http://cmtj-simulations.streamlit.app/)
![Downloads](https://img.shields.io/pypi/dm/cmtj.svg)

## Short description

A name may be misleading -- the MTJ (Magnetic Tunnel Junctions) are not the only structures that may be simulated.
The library allows for macromagnetic simulation of various multilayer spintronic structures. The package uses C++ implementation of (s)LLGS (stochastic Landau-Lifschitz-Gilbert-Slonczewski) equation with various field contributions included for instance: anisotropy, interlayer exchange coupling, demagnetisation, dipole fields etc.
It is also possible to connect devices in parallel or in series to have electrically coupled arrays.

## Demo

Check out the [streamlit hosted demo here](http://cmtj-simulations.streamlit.app/).

## Quickstart

#### Installation :rocket:

Installation is as easy as doing:
A recommended way is to use the `pip` package manager and virtualenv (or conda).

1. With `virtualenv` (recommended):

```bash
$(bash) python -m venv .my-venv
$(bash) source .my-venv/bin/activate
$(.my-venv) python -m pip install cmtj
```

2. Straight from `pip`:

```bash
python3 -m pip install cmtj
```

3. Straight from source:

```bash
python3 -m pip install https://github.com/LemurPwned/cmtj.git
```

4. Clone the repository:

```bash
git clone https://github.com/LemurPwned/cmtj.git
python3 -m pip install .
```

#### Extra dependencies

The package requires (if `utils` subpackage is used):

```
- numpy
- scipy
- matplotlib
```

## Read the docs

Documentation: [https://lemurpwned.github.io/cmtj](https://lemurpwned.github.io/cmtj)

## Extensions

There's a GUI version available! If you wish to conduct a subset of simulations, mainly for experimental modelling, please see the _PyMag_ project. It uses CMTJ as a backend for fast computation.

## Citing

We would appreciate citing either of the listed work if you decide to use the project or using the cite button on the right hand side panel of the repository:

[**cmtj: Simulation package for analysis of multilayer spintronic devices**](https://www.nature.com/articles/s41524-023-01002-x)

```bibtex
@article{mojsiejuk_cmtj_2023,
	title = {cmtj: Simulation package for analysis of multilayer spintronic devices},
	volume = {9},
	issn = {2057-3960},
	url = {https://www.nature.com/articles/s41524-023-01002-x},
	doi = {10.1038/s41524-023-01002-x},
	pages = {54},
	number = {1},
	journaltitle = {npj Comput Mater},
	author = {Mojsiejuk, Jakub and Ziętek, Sławomir and Grochot, Krzysztof and Skowroński, Witold and Stobiecki, Tomasz},
	date = {2023-04-06},
}
```

# Development

## Acknowledgements

Many thanks to professor Jack Sankey for his help with the development of thermal contributions, with inspiration from the [macrospinmob project](https://github.com/Spinmob/macrospinmob).

## Contributions

All contributions are welcome, please leave an issue if you've encountered any trouble with setup or running the library.

## Docker

In the `docker` directory there's a `Dockerfile` that can be used to build a docker image with the library installed.
`Dockerfile.app` is used for streamlit development.

## Precommit

There's a `.pre-commit-config.yaml` that does some basic python and cpp lints and checks. More static analysis to come in the future.
This may be run with:

```
pre-commit run -v
```

or

```
pre-commit run -a (or --files core/* cmtj/*)
```

## Documentation builds

There are couple of stages to building the documentation

1. Build Doxygen documentation
   ```
   doxygen Doxyfile
   ```
   This is mostly for the C++ documentation. Furture changes may couple C++ and Python docs.
2. Build stubs
   The stubgen is `pybind11-stubgen` or `mypy stubgen` with the latter being preferred now.
   E.g. to generate `Stack` module stubs we can go:
   ```
   stubgen -m cmtj.stack -o target-stub-dir/
   ```
   More info here: https://mypy.readthedocs.io/en/stable/stubgen.html.
3. Parse stubs to Markdown.
   This stage is done by running:
   `python3 docs/docgen.py `
   The deployment of the documentation is done via:
   ```bash
   mkdocs gh-deploy
   ```
   But first, worth a check:
   ```bash
   mkdocs serve
   ```
