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

- name: Khanak Bhargava
  affiliation: '1'
  
- name: Zhi Chen
  affiliation: '1'

- name: Eric Johnson
  affiliation: '1'

- name: Max P. Katz
  affiliation: '1'

- name: Piyush Sharda

- name:
  given-names: Alexander
  surname: Smith Clark
  affiliation: '1'

- name: Ben Wibking

- name: Donald Willcox

- name: Michael Zingale
  affiliation: '1'

affiliations:
- index: 1
  name: Department of Physics and Astronomy, Stony Brook University, Stony Brook, NY, USA

affiliations:

date: 01 January 2025

bibliography: paper.bib
---

# Summary

The AMReX-Astrophysics Microphysics library provides a common set of
microphysics routines (reaction networks, equations of state, and
other transport coefficients) for astrophysical simulation codes built
around the AMReX library.  Several codes, including Castro
[@castro_I], MAESTROeX [@maestroex], and Quokka [@quokka] use
Microphysics to provide the physics and solvers needed
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

Several classical Fortran libraries have been converted to header-only
C++ implementations, including the VODE integrator [@vode], the hybrid
Powell method of MINPACK [@powell], and the Runge-Kutta Chebyshev
integration method [@rkc].


# Design

Microphysics provides several different types of physics: equations of
state, reaction networks and screening methods, nuclear statistical
equilibrium solvers and table, conductivities, and opacities, as well
as the tools needed to work with them, most notably the ODE
integrators for the networks.

There are two ways to use Microphysics: standalone for simple investigations
or as part of an (AMReX-based) application code.  In both cases, the core
requirement is to select a network---this defines the composition that
is then used by most of the other physics routines.

A key design feature is the separation of the reaction network from
the integrator.  This allows us to easily experiment with different
integration methods (such as the RKC integrator) and also support
different modes of coupling reactions to a simulation code (operator
splitting and spectral deferred corrections [@castro_simple_sdc])

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

We also make use of the C++ autodiff library [@autodiff] to compute
thermodynamic derivatives required in the Jacobians of our reaction
networks.

# Capabilities


## equations of state

## networks

We have ported many of the classic "aprox" networks used in the
astrophysics community to C++ using templating to construct the
righthand side of the network at compile time.

Microphysics can also directly use networks created by the


# Unit tests / examples

Microphysics can be used as a standalone tool through the tests
in `Microphysics/unit_test/`.  There are 2 types of tests here:

* *comprehensive tests*: these test performance by setting up a cube
  of data (with density, temperature, and composition varying in a
  dimension) and performing an operation on the entire cube (calling
  the EOS, integrating a network, ...).

* *one-zone tests*: these simply call one of the physics modules with
  a single thermodynamic state.  This can be used to explore the
  physics that is implemented.

# References
