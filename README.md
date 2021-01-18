[![Build Status](https://travis-ci.org/starkiller-astro/Microphysics.svg?branch=main)](https://travis-ci.org/starkiller-astro/Microphysics) [![DOI](https://zenodo.org/badge/33425497.svg)](https://zenodo.org/badge/latestdoi/33425497)

# Microphysics

*A collection of astrophysical microphysics routines for stellar explosions*

There are several core types of microphysics routines hosted here:

* `EOS/`: these are the equations of state. All of them accept a struct
  called `eos_t` to pass the thermodynamic state information in
  and out, though in C++ they are templated such that they can accept
  other objects with members of the same name.

* `integration/`: this holds the various ODE integrators. VODE is the
  primary integrator for production use, but other integrators are provided
  for experimentation.

* `interfaces/`: this holds the structs used to interface with the
  EOS and networks.

* `networks/`: these are the reaction networks. They serve both to
  define the composition and its properties, as well as describe the
  reactions and energy release when reactions occur.

* `neutrinos/`: this holds the plasma neutrino cooling routines used
  in the reaction networks.

* `rates/`: this contains some common rate routines used by the
  various `aprox` networks, and could be expanded to contain other
  collections of rates in the future

* `screening/`: the screening routines for nuclear reactions. These
  are called by the various networks

* `unit_test/`: code specific to unit tests within this repo. In
  particular,

  - `test_eos` will test an equation of state by first calling
    it will (rho, T), and then calling it with other inputs
	to recover rho and/or T. A cube of data, with rho, T, and
	X is tested.

  - `test_react` will call a reaction network on a cube of
    data (rho, T, X).

  - `test_rhs` will evaluate the right-hand-side of a network
    without integrating it.

  - `test_sdc` tests integration using the set of variables for
    simplified SDC.

  - `test_jac` evaluates the Jacobian for comparison between the
    analytical and numerical Jacobians.

  - `test_screening` evaluates the screening routines.

* `util`: linear algebra routines for the various integrators
  (including BLAS and LINPACK), other math routines, and build
  scripts


# AMReX-Astro Codes

At the moment, these routines are written to be compatible with
the AMReX-Astro codes, Maestro and Castro.

* Castro: http://amrex-astro.github.io/Castro/

* Maestro: http://amrex-astro.github.io/MAESTRO/

To use this repository with AMReX codes, set `MICROPHYSICS_HOME` to
point to the `Microphysics/` directory.

There are various unit tests that work with the AMReX build system to
test these routines.


# Other Simulation Codes

The interfaces are fairly general, so they can be expanded to other
codes. This will require adding any necessary make stubs for the
code's build system as well as writing unit tests for that build
system to ensure the interfaces are tested.


# Documentation

A user's guide for Microphysics is available at:
http://starkiller-astro.github.io/Microphysics/docs/

The sphinx source for the documentation is in `Microphysics/sphinx_docs/`

## Development Model:

Development generally follows the following ideas:

  * New features are committed to the `development` branch.

    Nightly regression testing is used to ensure that no answers
    change (or if they do, that the changes were expected).

    If a change is critical, we can cherry-pick the commit from
    `development` to `main`.

  * Contributions are welcomed from anyone. *Any contributions that
    have the potential to change answers should be done via pull
    requests.*   A pull request should be generated from your fork of
    `Microphysics` and target the `development` branch. (If you mistakenly
    target `main`, we can change it for you.)

    Please add a line to `CHANGES` summarizing your change if it
    is a bug fix or new feature. Reference the PR or issue as
    appropriate. Additionally, if your change fixes a bug (or if
    you find a bug but do not fix it), and there is no current
    issue describing the bug, please file a separate issue describing
    the bug, regardless of how significant the bug is. If possible,
    in both the `CHANGES` file and the issue, please cite the pull
    request numbers or git commit hashes where the problem was
    introduced and fixed, respectively.

    If there are a number of small commits making up the PR, we may
    wish to squash commits upon merge to have a clean history.
    *Please ensure that your PR title and first post are descriptive,
    since these will be used for a squashed commit message.*

  * On the first workday of each month, we perform a merge of
    `development` into `main`, in coordination with `AMReX`,
    `Maestro`, and `Microphysics`. For this merge to take place, we
    need to be passing the regression tests.

    To accommodate this need, we close the merge window into
    `development` a few days before the merge day. While the merge
    window is closed, only bug fixes should be pushed into
    `development`. Once the merge from `development` -> `main` is
    done, the merge window reopens.


## Core Developers

People who make a number of substantive contributions will be named
"core developers" of StarKiller Microphysics. The criteria for
becoming a core developer are flexible, but generally involve one of
the following:

  * 10 non-merge commits to `Microphysics/` (including `Docs/`) or one
    of the problems that is not your own science problem *or*

  * addition of a new algorithm / module  *or*

  * substantial input into the code design process or testing

Core developers will be recognized in the following ways:

  * invited to the group's slack team

  * listed in the User's Guide and website as a core developer

  * invited to co-author general code papers / proceedings describing
    StarKiller Microphysics, its performance, etc. (Note: science
    papers will always be left to the science leads to determine
    authorship).

If a core developer is inactive for 3 years, we may reassess their
status as a core developer.
