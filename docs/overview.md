---
author:
  - LemurPwned
date: October 2022
title: Start here
---

# Overview

The `cmtj` library is a macrospin library. That means that it is efficient and easy to use. It is also a collection of utilities that are useful for macrospin analysis, such as selected procedures, energy calculations, and parallelization.

## Installation

Installation is as easy as doing:

```bash
python3 -m pip install cmtj
```

## How do I get started?

The best way to get started is to modify one of the examples in the [Examples](#CMTJBindingsTutorial) section.
If you're looking for a more in-depth look at the library, please see the [API](#API) section.

There are also some examples in the repository under [curated-examples folder](https://github.com/LemurPwned/cmtj/tree/master/curated-examples). These are less described

Remember to check out [#tips-and-trics](./tipsandtricks.md) section for some tips on how to improve your simulations.

### Functionalities

#### Junction system

Primarily you will be using the `Junction` class to create a multilayer device. This need not to be an MTJ, it can be of any ferromagnetic (FM) type. Remember, that the FM type is solely dependent on the parameters you pass to the constructor of the `Layer` class, such as $M_\textrm{s}$ or $K_\textrm{s}$ (magnetisation saturation and anisotropy respectively).

##### How do I model the spacer/oxide?

The funny thing is you don't do that directly. The effects of the spacer are effectively modeled by the `Junction` class in the following ways:

- the STT/SOT parameters, depending on the oxide type or thickness parameters such as $H_\textrm{DL}$ of $H_\textrm{FL}$ may get smaller or larger.
- the effective anisotropy in either of layers creating the junction may be dependent on the interface mixing or interface assymetry, so you can model that, for instance, by modifying one layer's anisotropy adequately.
- various modifications to magnetoresistance computation may be necessary, depending on the type of junction and the type of the oxide.
- addition of the `dipole` or `IEC` interaction is also a part of the junction model. Sometimes the coupling is so small, it may be neglected. See more in the [Core](#gen-docs/cmtj.md) section to see how to enable those (`setLayerIEC` and `setLayerTopDipole`).

#### Stack system

The junctions can be further stacked to create a stack of junctions. The stacking is done by creating a `Stack` object which can either be `Parallel` stack or `Series` stack. The type depends on the electrical connection between the junctions.

Find more in the examples.

#### Group interaction system

Another way to define 1-1 interaction matrix, also based on distance is `GroupInteraction` system. For example you can use it to define dipolar interaction in an array of devices. The interaction can also be custom, as long as it matches the signature -- for more, check out the [`interactions` tutorial](./tutorials/interactions.ipynb)

#### Symbolic models

There is a plethora of symbolic computation that can solve energy-based systems of FM layers, SAFs and so on.

For theory on how those models work visit [Smit-Beljers](./physics/sb-model.md) and [linearisation sections](./physics/linearisation.md). The corresponding API is given in [Symbolic models](./api/models/sb-general-reference.md) section.

Things they allow for:

- linearised frequency computation
- native Smit-Beljers model
- linearised Voltage Spin Diode

Check out the examples in [SBModel tutorial](./tutorials/SBModel.ipynb) and there are other examples scattered here and there.
