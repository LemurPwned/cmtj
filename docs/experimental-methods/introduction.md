---
author:
  - LemurPwned
date: April 2024
title: Introduction to experimental methods
---

# Experiment modelling

`cmtj` is a versatile tool for simulating a wide range of experiments.
It's possible to simulate a variety of experiments, from simple to complex, below is the list of the things we've tried so far and verified:

- R(H) and M(H) loops
- FMR -- ferromagnetic resonance, both field, current and frequency sweeps.
  - SD-FMR -- spin diode FMR (VSD)
  - STT-FMR -- spin transfer torque FMR
  - SOT-FMR -- spin orbit torque FMR
- PIMM -- pulse induced magnetometry
- Harmonic Hall voltage (HHV) measurements
  - Field magnitude HHV (at a set angle)
  - Angular HHV (a function of field angle)
- CIMS -- current induced magnetisation switching

For those experiments we have many examples, along with some data to compare to:

- [R(H) and M(H) loops](../examples/#mh-and-rh-loops)
- [Angular Harmonic Hall voltage](../examples/#angular-harmonics) and [Field magnitude Harmonic Hall voltage](../HarmonicsFits)
- [CIMS](../examples/#cims) and [here](../CIMS)
- Spin Diode [FMR](../examples/#spin-diode-fmr) and [here](../VoltageSpinDiodeFits)

!!! note

      These are the ones we use in our lab, but that does not mean that you can't simulate other experiments. Let us know if you have a specific experiment in mind to replicate and fit and we'll try to help you out.
