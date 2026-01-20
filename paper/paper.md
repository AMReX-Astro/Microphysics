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
- surname: AMReX-Astro Microphysics Team
  affiliation: '†'

- name: Khanak Bhargava
  affiliation: '1'

- name: Abigail Bishop
  affiliation: '2'

- name: Zhi Chen
  affiliation: '1'

- name: Doreen Fan
  affiliation: '3'

- name: Carl Fields
  affiliation: '4'

- name: Adam M. Jacobs
  affiliation: '3'
  
- name: Eric T. Johnson
  affiliation: '1'

- name: Max P. Katz
  affiliation: '1'

- name: Mark Krumholz
  affiliation: '5'
  
- name: Chris Malone
  affiliation: '6'

- name: Andy Nonaka
  affiliation: '7'
  
- name: Piyush Sharda
  affiliation: '8'

- given-names: Alexander
  surname: Smith Clark
  affiliation: '1'

- name: Frank Timmes
  affiliation: '9'

- name: Ben Wibking
  affiliation: '10'

- name: Don E. Willcox
  affiliation: '3'
  
- name: Michael Zingale
  affiliation: '1'

affiliations:
- index: †
  name: https://github.com/amrex-astro/Microphysics
- index: 1
  name: Department of Physics and Astronomy, Stony Brook University, Stony Brook, NY, USA
- index: 2
  name: Department of Physics, University of Wisconsin, Madison, Madison, WI, USA
- index: 3
  name: affiliation not disclosed
- index: 4
  name: Department of Astronomy, University of Arizona, Tucson, AZ, USA
- index: 5
  name: Research School of Astronomy and Astrophysics, The Australian National University, Australia
- index: 6
  name: Los Alamos National Laboratory, Los Alamos, NM, USA
- index: 7
  name: Lawrence Berkeley National Laboratory, Berkeley, CA, USA
- index: 8
  name: Leiden Observatory, Leiden, The Netherlands
- index: 9
  name: Arizona State University, Tempe, AZ, USA
- index: 10
  name: Department of Physics and Astronomy, Michigan State University, E. Lansing, MI, USA

date: 20 July 2025

bibliography: paper.bib
---


# Summary

The AMReX-Astrophysics Microphysics library provides a common set of
microphysics routines (reaction networks and associated physics,
equations of state, and various transport coefficients) as well as
solvers (stiff ODE integrators, nonlinear system solvers) for
astrophysical simulation codes built around the AMReX adaptive mesh
refinement library [@amrex].  Several multi-dimensional simulation
codes, including the compressible hydrodynamics code Castro
[@castro_I], the low-Mach number hydrodynamics code MAESTROeX
[@maestroex], and the radiation-hydrodynamics code Quokka [@quokka]
use Microphysics to provide the physics and solvers needed to close
the hydrodynamics systems that they evolve.  The library is
implemented in C++ with GPU-offloading a key design feature.

# Statement of need

Astrophysical simulation codes need many different smallscale
(microphysics) physics inputs to close the system of equations.  There
are many astrophysics simulation codes built around the AMReX library,
with each specializing in different astrophysics phenomena.  Each of
these codes share some common needs.  The Microphysics library was
created to minimize developer effort across these codes and coordinate
the approach to exascale compute architectures, in particular, GPU
support for astrophysical simulation codes.


# Design

The Microphysics project started in 2013 as a way to centralize the
reaction networks and equations of state used by Castro and MAESTRO
[@maestro], the predecessor to MAESTROeX.  Originally, Microphysics
used Fortran and for a brief period, it was referred to as Starkiller
Microphysics, which was an attempt to co-develop microphysics routines
for the Castro and the Flash [@flash] simulation codes.  As interest
in GPUs grew (with early support added to Microphysics in 2015),
Castro moved from a mix of C++ and Fortran to pure C++ to take
advantage of GPU-offloading afforded by the AMReX library, and C++
ports of all physics routines and solvers were added to Microphysics.
At this point, the project was formally named the AMReX-Astrophysics
Microphysics library.  Today, the library is completely written in C++
and relies heavily on the AMReX data structures to take advantage of
GPUs.  The GPU-enabled reaction network integrators led to the Quokka
code adopting Microphysics for their simulations.

Microphysics provides several different types of physics: equations of
state, reaction networks and screening methods, nuclear statistical
equilibrium solvers and tabulations, thermal conductivities, and
opacities, as well as the tools needed to work with them, most notably
the suite of stiff ODE integrators for the networks.
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

Microphysics uses header-only implementations of all functionality as
much as possible, to allow for easier compiler inlining, which is
especially important in GPU kernels.  We also leverage C++17 `if
constexpr` templating to compile out unnecessary computations for
performance.  Generally, the physics routines and solvers are written
to work on a single zone from a simulation code, and in AMReX, a C++
lambda-capturing approach is used to loop over zones (and offload to
GPUs if desired).  When used with an application code, this design
permits the simulation state data to be allocated directly in GPU
memory and left there for the entire simulation, with all physics run
directly on the GPU.  Since each zone in a simulation usually will
have a different thermodynamic state, the integration of reaction
networks can lead to thread divergence issues, since some zones will
have an easier burn than others.  To help mitigate this issue, we can
cap the number of integration steps and either retry an integration on
a zone-by-zone basis with different tolerances or Jacobian
approximations or pass the failure back to the application code to
deal with.  This strategy has been successful for many large scale
simulations [@Zingale_2025].


Another key design feature is the separation of the reaction network
from the integrator.  This allows us to easily experiment with
different integration methods (such as the RKC integrator) and also
support different modes of coupling reactions to a simulation code,
including operator splitting and spectral deferred corrections (SDC)
(see, e.g., @castro_simple_sdc).  The latter is especially important
for explosive astrophysical flows.  Tight integration with pynucastro [@pynucastro], [@pynucastro2], allows for the generation of custom reaction networks for a science problem.

There are two ways to use Microphysics: in a standalone fashion (via
the unit tests) for simple investigations or as part of an
(AMReX-based) application code.  In both cases, the core
(compile-time) requirement is to select a network---this defines the
composition that is then used by most of the other physics routines.
This choice is done at compile-time, allowing
Microphysics to provide the number of species as a `constexpr` value
(which many application codes need), and also greatly reducing the
compilation time (due to the templating used throughout the library).

# Research Impact Statement

Microphysics has been used for simulations of convective Urca
[@Boyd_2025] and X-ray bursts [@Guichandut_2024] with MAESTROeX; and
for simulations of nova [@Smith2025], X-ray bursts [@Harpole_2021],
thermonuclear supernovae [@Zingale_2024_dd], and convection in massive
stars [@Zingale_2024] with Castro. This Microphysics library has also
enabled recent work in astrophysical machine learning to train deep
neural networks modeling nuclear reactions [@nn_astro_2022], [@dnn_astro_2025].

# AI Usage Disclosure

No generative AI/LLM was used for code or documentation generation in the
git repository.  We have experimented with using AI/LLM tools for code
review and for suggesting places to focus our optimization efforts on,
but the resulting coding, benchmarking, and testing is then done by
humans.


# Acknowledgements

The AMReX-Astro Microphysics library developers are an open scientific
team with members contributing to various aspects of the library. We
have thus chosen to display the members of our development team in the
author list in alphabetical order.  All developers who have
contributed new features, substantial design input, and/or at least 3
commits were invited to be coauthors.  The work at Stony Brook was
supported by the US Department of Energy, Office of Nuclear Physics
grant DE-FG02-87ER40317.

# References

