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
  affiliation:
  
- name: Khanak Bhargava
  affiliation: '1'
  
- name: Zhi Chen
  affiliation: '1'

- name: Eric Johnson
  affiliation: '1'

- name: Max P. Katz
  affiliation: '1'

- name: Piyush Sharda
  affiliation: '2'
  
- given-names: Alexander
  surname: Smith Clark
  affiliation: '1'

- name: Ben Wibking
  affiliation: '3'

- name: Donald Willcox

- name: Michael Zingale
  affiliation: '1'

affiliations:
- index: 1
  name: Department of Physics and Astronomy, Stony Brook University, Stony Brook, NY, USA
- index: 2
  name: Leiden Observatory, Leiden, The Netherlands
- index: 3
  name: Department of Physics and Astronomy, Michigan State University, E. Lansing, MI, USA

date: 01 January 2025

bibliography: paper.bib
---


# Summary

The AMReX-Astrophysics Microphysics library provides a common set of
microphysics routines (reaction networks, equations of state, and
other transport coefficients) for astrophysical simulation codes built
around the AMReX adaptive mesh refinement library [@amrex].  Several
multi-dimensional simulation codes, including the compressible hydrodynamics code Castro
[@castro_I], the low-Mach number hydrodynamics code MAESTROeX
[@maestroex], and the radiation-hydrodynamics code Quokka [@quokka]
use Microphysics to provide the physics and solvers needed to close
the hydrodynamics systems that they evolve.

# History

This project in started out in 2013 as a way to centralize the
reaction networks and equations of state used by Castro and MAESTRO
[@maestro], the predecessor to MAESTROeX.  Originally, Microphysics used Fortran and 
for a brief period, it was
referred to as Starkiller Microphysics, which was an attempt to
co-develop microphysics routines for the Castro and the Flash [@flash]
simulation codes.   As interest in GPUs grew (with early support added to Microphysics in 2015),
Castro moved from a mixed of C++ and Fortran to pure C++ to take
advantage of GPU-offloading afforded by the AMReX library and C++ ports were
added to Microphysics.  At this point,
the development focused solely on AMReX-based codes and C++ and the project
was formally named the AMReX-Astrophysics Microphysics library and the
Fortran implementations were removed over time.
Today, the library is completely written in C++ and relies heavily on
the AMReX data structures to take advantage of GPUs.

# Design

Microphysics provides several different types of physics: equations of
state, reaction networks and screening methods, nuclear statistical
equilibrium solvers and tabulations, thermal conductivities, and
opacities, as well as the tools needed to work with them, most notably
the suite of stiff ODE integrators for the networks.

There are two ways to use Microphysics: standalone for simple investigations
or as part of an (AMReX-based) application code.  In both cases, the core
requirement is to select a network---this defines the composition that
is then used by most of the other physics routines.

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

Several classic Fortran libraries have been converted to header-only
C++ implementations, including the VODE integrator [@vode], the hybrid
Powell method of MINPACK [@powell], and the Runge-Kutta Chebyshev
(RKC) integration method [@rkc].  The code was modernized where possible,
with many `go to` statements removed and additional logic added
to support our applications (see for example the discussion
on VODE in @castro_simple_sdc).
We also make use of the C++ autodiff library [@autodiff] to compute
thermodynamic derivatives required in the Jacobians of our reaction
networks.

Another key design feature is the separation of the reaction network
from the integrator.  This allows us to easily experiment with
different integration methods (such as the RKC integrator) and also
support different modes of coupling reactions to a simulation code,
including operator splitting and spectral deferred corrections (SDC)
(see, e.g., @castro_simple_sdc).  The latter is especially important
for explosive astrophysical flows.


# Capabilities


## Reaction networks

A reaction network defines the composition (including the atomic
weight and number) and the reactions that link the nuclei together.
Even if reactions are not being modeled, a `general_null` network can
be used to simply define the composition.

In multidimensional simulations, there is a desire to make the
reaction as small as possible (due to the memory and per-zone
computational costs) while still being able to represent the
nucleosynthesis reasonable accurately.  As a result, approximations
to rates are common and a wide variety of networks are used depending
on the burning state being modeled.

We have ported many of the classic "aprox" networks used in the
astrophysics community (for example "aprox21" described in
@wallacewoosley:1981 to C++.  Many of these originated from the
implementations of @cococubed.  Our implementation relies heavily on
C++ templates, allowing us to simply define the properties of the
reactions and then the compiler builds the righthand side and Jacobian
of the system at compile-time.  This reduces the maintenance costs of
the networks and also eliminates some common indexing bugs.

We also integrate with the pynucastro nuclear astrophysics library
[@pynucastro; @pynucastro2], allowing us to generate a custom network
in a few lines of python simply by specifying the nuclei we want.  This
makes use of the reaction rates from @ReacLib and others, and allows us
to keep up to date with changes in rates and build more complex networks
than the traditional aprox nets.


## Equations of state




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
