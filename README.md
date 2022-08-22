<img style="float: right; max-width: 50px;" src="docs/assets/icon.svg">

# CMTJ

[![PyPI](https://github.com/LemurPwned/cmtj/actions/workflows/main.yml/badge.svg?branch=master)](https://github.com/LemurPwned/cmtj/actions/workflows/main.yml)
[![pages-build-deployment](https://github.com/LemurPwned/cmtj/actions/workflows/pages/pages-build-deployment/badge.svg?branch=gh-pages)](https://github.com/LemurPwned/cmtj/actions/workflows/pages/pages-build-deployment)
[![Version](https://img.shields.io/pypi/v/cmtj)]()
[![License](https://img.shields.io/pypi/l/cmtj.svg)]()
[![Downloads](https://img.shields.io/pypi/dm/cmtj.svg)]()

## Short description

The library allows for macromagnetic simulation of multilayer spintronic structures.
A name may be misleading -- the MTJ (Magnetic Tunnel Junctions) are not the only structures that may be simulated.

## Quickstart

#### Installation :rocket:

Installation is as easy as doing:
A recommended way is to use the `pip` package manager and virtualenv (or conda).

1. With `virtualenv`

```bash
$(bash) python3 -m venv .my-venv
$(bash) source .my-venv/bin/activate
$(.my-venv) python -m install cmtj
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

#### Read the docs

Documentation: [https://lemurpwned.github.io/cmtj](https://lemurpwned.github.io/cmtj)

## WIKI :mortar_board:

Read more in do [the docs here](https://lemurpwned.github.io/cmtj/).

## Extensions

There's a GUI version available! If you wish to conduct a subset of simulations, mainly for experimental modelling, please see the _PyMag_ project. It uses CMTJ as a backend for fast computation.

## Citing

Please cite if you decide to use the project

```bibtex
@article{PhysRevB.106.024403,
  title = {Numerical model of harmonic Hall voltage detection for spintronic devices},
  author = {Zi\ifmmode \mbox{\k{e}}\else \k{e}\fi{}tek, S\l{}awomir and Mojsiejuk, Jakub and Grochot, Krzysztof and \L{}azarski, Stanis\l{}aw and Skowro\ifmmode \acute{n}\else \'{n}\fi{}ski, Witold and Stobiecki, Tomasz},
  journal = {Phys. Rev. B},
  volume = {106},
  issue = {2},
  pages = {024403},
  numpages = {9},
  year = {2022},
  month = {Jul},
  publisher = {American Physical Society},
  doi = {10.1103/PhysRevB.106.024403},
  url = {https://link.aps.org/doi/10.1103/PhysRevB.106.024403}
}
```

# Development

## Acknowledgements

Many thanks to professor Jack Sankey for his help with the development of thermal contributions, with inspiration from the [macrospinmob project](https://github.com/Spinmob/macrospinmob).

## Contributions

All contributions are welcome, please leave an issue if you've encountered any trouble with setup or running the library.

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
