# CMTJ Curated Examples

This directory contains high-quality, well-documented examples demonstrating key capabilities of the CMTJ simulation package for spintronic devices.

## Overview

These examples follow best practices from `AGENTS.md` and serve as the golden standard for CMTJ simulations. Each example includes:
- Comprehensive documentation explaining the physics and mechanisms
- Proper parameter values following recommended ranges
- Correct unit handling (Ms in T for core layers, A/m for domain walls/SB models)
- Complete simulation workflow from setup to visualization

## Example Categories

### Magnetization Dynamics (LLG-based)

1. **switching-hysteresis.py** - Basic magnetization switching and hysteresis loops
   - Single layer field-driven switching
   - Hysteresis behavior and coercive field extraction
   - Good starting point for beginners

2. **sto.py** - Spin-Transfer Oscillator
   - Persistent magnetization oscillations via spin-transfer torque
   - Trajectory visualization on unit sphere
   - STT dynamics demonstration

3. **p-bit.py** - Probabilistic Bit with Thermal Noise
   - Stochastic magnetization switching
   - Statistical analysis of switching events
   - Thermal activation and autocorrelation

4. **fmr-basic.py** - Ferromagnetic Resonance Spectroscopy
   - Multi-layer FMR with pulsed excitation
   - Frequency-field spectrograms
   - FFT analysis of magnetization dynamics

### Coupling and Interactions

5. **iec-coupling.py** - Interlayer Exchange Coupling
   - Two-layer system with IEC interaction
   - Antiferromagnetic coupling demonstration
   - Coupled precession modes

6. **dw-motion.py** - Domain Wall Motion
   - Current-driven domain wall dynamics
   - SOT effects and DMI stabilization
   - Collective coordinate model

### Current-Induced Magnetization Switching (CIMS)

7. **cims-hysteresis.py** - CIMS Hysteresis Loops
   - Spin-orbit torque switching
   - Current-field phase diagrams

8. **cims-stability-stability-diagram.py** - CIMS Stability Diagrams
   - Comprehensive SOT stability analysis
   - Parameter optimization for memory devices

### Experimental Methods

9. **harmonic-hall.py** - Harmonic Hall Voltage Analysis
   - Second harmonic Hall measurements
   - SOT characterization technique

10. **ahe_loops.py** - Anomalous Hall Effect Loops
    - Resistance vs field measurements
    - Layer-resolved switching detection

11. **vsd-basic-field-scan.py** - Voltage Spin Diode (Field Scan)
    - DC voltage generation from RF currents
    - Field-dependent mixing effects

12. **vsd-basic-freq-scan.py** - Voltage Spin Diode (Frequency Scan)
    - Frequency-resolved spin diode response
    - Resonance detection

## Quick Start

1. Install CMTJ with dependencies:
```bash
pip install -e .
pip install matplotlib scipy scienceplots
```

2. Run any example:
```bash
python curated-examples/switching-hysteresis.py
```

3. Check the `figures/` subdirectory for output plots.

## Parameter Guidelines

All examples follow the parameter ranges from `AGENTS.md`:

| Parameter | Typical Range | Units | Notes |
|-----------|---------------|-------|-------|
| Ms | 0.5 – 1.6 | T | Use A/m for DW/SB models |
| Ku | 1e2 – 1e6 | J/m³ | PMA range |
| α (damping) | 0.01 – 0.03 | - | Gilbert damping |
| J (IEC) | ±0.001 – ±3.0 | mJ/m² | Interlayer coupling |
| D (DMI) | 0 – 3.0 | mJ/m² | Interfacial DMI |
| H_ext | ±0 – ±500e3 | A/m | External field |
| dt | 1e-12 | s | Time step (unless adaptive) |
| sim_time | 1e-9 – 500e-9 | s | Simulation duration |

## New Examples Added

### 1. Domain Wall Motion (dw-motion.py)
Demonstrates domain wall dynamics in magnetic nanowires:
- Spin-orbit torque driving
- DMI stabilization of Néel walls
- Position, velocity, and width evolution
- Pinning effects

**Key Physics:** Collective coordinate approach to DW motion, SOT-driven dynamics

### 2. IEC Coupling (iec-coupling.py)
Shows interlayer exchange coupling between magnetic layers:
- Antiferromagnetic coupling (J1 < 0)
- Coupled precession modes
- Phase space analysis
- Field-driven dynamics

**Key Physics:** RKKY-type coupling, acoustic/optical modes

### 3. Switching and Hysteresis (switching-hysteresis.py)
Basic field-driven magnetization reversal:
- Hysteresis loop generation
- Field sweep protocols
- Coercive field extraction
- Magnetization trajectories

**Key Physics:** Energy barrier crossing, hysteretic behavior

## Contributing

When adding new examples:
1. Follow the documentation style of existing examples
2. Include comprehensive physics explanation in docstring
3. Use recommended parameter ranges from AGENTS.md
4. Generate and save a figure showing results
5. Test that the example runs successfully
6. Add entry to this README

## References

- Main documentation: `docs/`
- Physics background: `docs/physics/contributions.md`
- API reference: `docs/api/`
- Agent guidelines: `AGENTS.md`
