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

- name: Piyush Sharda
  affiliation: '7'

- given-names: Alexander
  surname: Smith Clark
  affiliation: '1'

- name: Frank Timmes
  affiliation: '8'

- name: Ben Wibking
  affiliation: '9'

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
  name: Leiden Observatory, Leiden, The Netherlands
- index: 8
  name: Arizon State University, Tempe, AZ, USA
- index: 9
  name: Department of Physics and Astronomy, Michigan State University, E. Lansing, MI, USA

date: 13 July 2025

bibliography: paper.bib
---


# Summary

The AMReX-Astrophysics Microphysics library[^authors] provides a common set of
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
the approach to exascale compute architectures.  Microphysics has been
used for simulations of convective Urca [@Boyd_2025] and X-ray bursts
[@Guichandut_2024] with MAESTROeX; and for simulations of nova
[@Smith2025], X-ray bursts [@Harpole_2021], thermonuclear supernovae
[@Zingale_2024_dd], and convection in massive stars [@Zingale_2024]
with Castro. This Microphysics library has also enabled recent work
in astrophysical machine learning to train deep neural networks
modeling nuclear reactions in [@nn_astro_2022] and [@dnn_astro_2025].

# Project history

The Microphysics project started in 2013 as a way to centralize the
reaction networks and equations of state used by Castro and MAESTRO
[@maestro], the predecessor to MAESTROeX.  Originally, Microphysics
used Fortran and for a brief period, it was referred to as Starkiller
Microphysics, which was an attempt to co-develop microphysics routines
for the Castro and the Flash [@flash] simulation codes.  As interest
in GPUs grew (with early support added to Microphysics in 2015),
Castro moved from a mix of C++ and Fortran to pure C++ to take
advantage of GPU-offloading afforded by the AMReX library and C++
ports of all physics routines and solvers were added to Microphysics.
At this point, the development focused solely on AMReX-based codes and
C++ and the project was formally named the AMReX-Astrophysics
Microphysics library.  Today, the library is completely written in C++
and relies heavily on the AMReX data structures to take advantage of
GPUs.  The GPU-enabled reaction network integrators led to the Quokka
code adopting Microphysics for their simulations.

# Design

Microphysics provides several different types of physics: equations of
state, reaction networks and screening methods, nuclear statistical
equilibrium solvers and tabulations, thermal conductivities, and
opacities, as well as the tools needed to work with them, most notably
the suite of stiff ODE integrators for the networks.

There are two ways to use Microphysics: in a standalone fashion (via
the unit tests) for simple investigations or as part of an
(AMReX-based) application code.  In both cases, the core
(compile-time) requirement is to select a network---this defines the
composition that is then used by most of the other physics routines.

Microphysics uses header-only implementations of all functionality as
much as possible, to allow for easier compiler inlining.  Generally,
the physics routines and solvers are written to work on a single zone
from a simulation code, and in AMReX, a C++ lambda-capturing approach
is used to loop over zones (and offload to GPUs if desired).  We also
leverage C++17 `if constexpr` templating to compile out unnecessary
computations for performance.  For example, our equations of state can
compute a lot of thermodynamic quantities and derivatives, but for
some operations, we only need a few of these.  All of the equations of
state are templated on the `struct` that holds the thermodynamic
state.  If we pass the general `eos_t` type into the EOS, then
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

Finally, most of the physics is chosen at compile-time.  This allows
Microphysics to provide the number of species as a `constexpr` value
(which many application codes need), and also greatly reduces the
compilation time (due to the templating used throughout the library).

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
@wallacewoosley:1981) to C++.  Many of these originated from the
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


### Screening

Nuclear reaction rates are screened by the electrons in the plasma
(which reduce the Coulomb barrier for the positively charged nuclei to
fuse).  Microphysics provides several different screening
implementations: the widely-used `screen5` method based on
[@graboske:1973; @jancovici:1977; @alastuey:1978; @itoh:1979], the
methods of [@chugunov:2007] and [@chugunov:2009], and the method of
[@Chabrier_1998].


### Nuclear statistical equilibrium

At high temperatures ($T > 4\times 10^9~\mathrm{K}$), forward and
reverse reactions can come into equilibrium (nuclear statistical
equilibrium, NSE).  Integrating the reaction network directly in this
regime can be difficult, since the large, but oppositely signed rates,
may not cancel exactly.  In this case, instead of integrating the
network, we can impose the equilibrium state.  Microphysics has two
different approaches to NSE: a self-consistent solve for the NSE state
using the nuclei in the present reaction network (similar to
@Kushnir_2020) and an interpolation from a tabulated NSE state that
was generated with $\mathcal{O}(100)$ nuclei (see @Zingale_2024).

### Thermal neutrinos

There are a number of thermal mechanisms for producing neutrinos,
including plasma, photo, pair, recombination, and Bremsstrahlung
neutrinos.  These act as an energy loss term to the reaction network
and are implemented following @itoh:1996.




## Equations of state

The equations of hydrodynamics are closed via an equation of state
that related internal energy, pressure, and density (along with
composition).  For systems with reactions or thermal diffusion, it
also provides temperature.  Traditionally, equations of state are
implemented in terms of density and temperature, so a Newton-Raphson
method is used to invert the EOS given energy and density (or some
other thermodynamic quantities).  A wide range of thermodynamic
quantities are needed by simulation codes, including pressure,
internal energy, enthalpy, entropy, and their derivatives with
respect to density, temperature, and composition.  The various EOS
`struct` types carry this thermodynamic state.

A variety of EOSs are implemented, to allow for application to a range
of problems.  These include a simple gamma-law EOS, the stellar EOS of
@timmes:2000, and an equation of state applicable to primordial
chemistry.

## Transport coefficients

For thermal diffusion or radiation transport, conductivities and
opacities are needed.  We provide a C++ port of the stellar
conductivity opacities from @timmes:2000b.  These are appropriate for
modeling thermonuclear flames in supernovae and X-ray bursts.

# GPU Strategy

Microphysics is designed such that all computation takes place on
GPUs.  When used with an application code, this permits the simulation
state data to be allocated directly in GPU memory and left there for
the entire simulation.  For the ODE integration, the integrator
(e.g. VODE) is run on the GPU directly.  Since each zone in a
simulation usually will have a different thermodynamic state, this can
lead to thread divergence issues, since some zones will have an easier
burn than others.  To help mitigate this issue, we can cap the number
of integration steps and either retry an integration on a zone-by-zone
basis with different tolerances or Jacobian approximations or pass the
failure back to the application code to deal with.  This strategy
has been successful for many large scale simulations [@Zingale_2025].


# Unit tests / examples

Microphysics can be used as a standalone tool through the tests
in `Microphysics/unit_test/`.  There are 2 types of tests here:

* *comprehensive tests*: these test performance by setting up a cube
  of data (with density, temperature, and composition varying in a
  dimension) and performing an operation on the entire cube (calling
  the EOS, integrating a network, ...).  A separate test is provided
  for each major physics module.

* *one-zone tests*: these simply call one of the physics modules with
  a single thermodynamic state.  This can be used to explore the
  physics that is implemented, and also serve to demonstrate the interfaces
  used in Microphysics.

These tests also serve as tutorial codes for integrating Microphysics
into new application codes.

# Acknowledgements

All developers who have contributed new features, substantial design
input, and/or at least 3 commits were invited to be coauthors.  The
work at Stony Brook was supported by the US Department of Energy,
Office of Nuclear Physics grant DE-FG02-87ER40317.

# References

[^authors]: The AMReX-Astro Microphysics library developers are an open scientific team with members contributing to various aspects of the library. We have thus chosen to display the members of our development team in the author list in alphabetical order.
