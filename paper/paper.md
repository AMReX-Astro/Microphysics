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

bibliography: refs.bib
---

# Summary

The AMReX-Astrophysics Microphysics library provides a common set of
microphysics routines (reaction networks, equations of state, and
other transport coefficients) for astrophysical simulation codes built
around the AMReX library.  Several codes, including Castro, MAESTROeX,
and Quokka use Microphysics to provide the physics and solvers needed
to close the hydrodynamics systems that they evolve.  

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

Several classical Fortran libraries have been converted to header-only C++
implementations, including the VODE integrator and the hybrid Powell method
of MINPACK.  


# Design

Microphysics provides several different types of physics: equations of
state, reaction networks and screening methods, nuclear statistical
equilibrium solvers and table, conductivities, and opacities, as well
as the tools needed to work with them, most notably the ODE integrators
for the networks.

There are two ways to use Microphysics: standalone for simple investigations
or as part of an (AMReX-based) application code.  In both cases, the core
requirement is to select a network---this defines the composition that
is then used by most of the other physics routines.

A key design feature is the separation of the reaction network from
the integrator.  This allows us to easily experiment with different
integration methods (such as the RKC integrator) and also support
different modes of coupling reactions to a simulation code (operator
splitting and spectral deferred corrections)

We rely on header-only implementations as much as possible, to allow
for easier compiler inlining.  We also leverage C++17 `if constexpr`
templating to compile out unnecessary computations for performance.
For example, our equations of state can compute a lot of thermodynamic
quantities and derivatives, but for some operations, we only need a
few of these.  If we pass the general `eos_t` type into the EOS, then
everything is calculated, but if we pass in to the same interface the
smaller `eos_re_t` type, then only a few energy terms are computed
(those that are needed when finding temperature from specific internal
energy).


pynucastro integration


# Capabilities




# Unit tests / examples

Microphysics can be used as a standalone tool through the tests
in `Microphysics/unit_test/`.  There are 2 types of tests here:

* *comprehensive tests* test performance by setting up a cube of data
  (with density, temperature, and composition varying in a dimension)
  and performing an operation on the entire cube (calling the EOS,
  integrating a network, ...).

* *one-zone tests* simply call one of the physics modules with a
  single thermodynamic state.  This can be used to explore the physics
  that is implemented.

# References

