---
title: 'AMReX-Astrophysics Microphysics: A set of microphysics routines for astrophysical simulation codes
based on the AMReX library'

tags:
- nuclear
- reactions
- astrophysics
- physics
- integration
- differential equations

authors:
- name: AMReX-Astro Microphysics Development Team

affiliations:

date: 01 January 2025

bibliography: paper.bib
---

# Summary

The AMReX-Astrophysics Microphysics library provides a common set of
microphysics routines (reaction networks, equations of state, and
other transport coefficients) for astrophysical simulation codes built
around the AMReX library.

# History

This project started out as Starkiller Microphysics, which was an
attempt to codevelop microphysics routines for the Castro and Flash
simulation codes.  Originally the library used Fortran and was
restricted to CPUs, but C++ ports were added over time to take
advantage of GPU-offloading afforded by the AMReX library.  Eventually
as the development of the two codes diverged, the C++ ports of the
Microphysis were split off into the AMReX-Astrophysics Microphysics
library.  Today, the library is completely written in C++ and relies
on the AMReX data structures.  

# Design

Porting many classic solvers from old Fortran to C++

header-only library as much as possible

C++17

pynucastro integration


# Capabilities



# Unit tests / examples


# References

