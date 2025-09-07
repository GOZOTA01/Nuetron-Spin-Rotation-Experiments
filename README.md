# Neutron Spin Rotation Experiments

A Fortran-based simulation suite for neutron spin transport through magnetic field systems, designed for magnetic spectrometer and neutron beam transport analysis.

## Overview

This project simulates neutron trajectories and spin precession through complex magnetic field geometries using quantum mechanical rotation operators and Larmor precession physics.

## Features

- **Three Simulation Modes:**
  - Mode 0: Quadrant Simulation (4-point symmetry test)
  - Mode 1: Numsteps Simulation (systematic spatial grid scan)
  - Mode 2: Random Simulation (Monte Carlo statistical analysis)

- **Advanced Physics:**
  - Neutron spin transport through magnetic fields
  - Larmor precession calculations
  - Adiabatic spin following
  - Quantum mechanical rotation operators
  - 3D magnetic field interpolation

## Files

- `Outcoil.f` - Main program controlling simulation workflow
- `sBrot.f` - Core neutron spin transport calculations
- `Bfield.f` & `Bfield2.f` - Magnetic field calculation routines
- `trilinearinterp.f` - 3D field interpolation
- `makefile` - Build system
- `heatmap.py` - Python visualization script

## Building

```bash
make clean
make
