changes since last release:

   * VODE90 now works with the simplified SDC time step algorithms,
     and the preprocessor option for this SDC was changed to
     SIMPLIFIED_SDC (#194)

# 19.10

   * The iso7 network was ported to GPUs (#172)

   * VODE90 now better agrees with VODE (#192)

   * When building with multiple integrators, the contents of the rpar
     modules could clash.  This has been fixedc (#136)

   * A module for making the "Nonaka plot" tracking the evolution of a
     quantity during the burn was added, and is enabled with
     USE_NONAKA_PLOT=TRUE

# 19.08

   * A new network, subch2, was added that combines aprox13 and the
     subch networks. (#184)

# 19.05

   * The aprox21 network was missing the analytic Jacobian term for
     the derivative of He4 with respect to Ni56. This is fixed. (#175)

   * The numerical Jacobian module, used by the BS and VBDF integrators
     had some wrong scalings.  These have now been fixed (#179, #180)

# 19.01

  * the docs are now automatically build from the sphinx source
    using travis on github.

# 18.12

  * simplify conductivity interface to match the eos interface by
    moving the conductivity into the eos type

# 18.11

  * new python bindings have been added

  * the documentation has been switched to Sphinx and is now hosted
    online.

  * a bug was fixed in the stellarcollapse EOS where the energy offset
    was being applied incorrectly.


# 18.10

  * test_eos and test_react now both work on GPUs (using the AMReX
    `gpu` branch)

  * the intermediate blending of the weak and strong screening regimes
    was wrong, and has been fixed.  We've also synced some parameters
    up to agree with those in MESA and Kepler.  (#149, 150)

  * eos_input_is_constant is now set to true for the helmholtz EOS.
    This mean that the EOS inputs will not be modified after the EOS
    call.  This is good for conserving energy in a hydro code, but the
    tradeoff is a small (to root finding tolerance) inconsistency in
    the thermodynamic state.

# 18.09

  * The Helmholtz parameters ttol and dtol (controlling the error
    for the Newton iteration when in a mode other than eos_input_rt)
    are now runtime parameters in the extern namelist as eos_ttol
    and eos_dtol.

# 18.08

  * the unit tests (test_react, test_sdc, and test_eos) have been
    ported from the Fortran to C++ build system in AMReX.  This will
    allow us to test the GPU framework in AMReX.

# 18.07

  * added CUDA support to the VODE90 integrator, the helmeos, and the
    networks aprox13, aprox19, aprox21, ignition_simple,
    C-burn-simple, URCA-simple.

  * Ported the unit test frameworks to FBoxLib

# 18.05

  * lots of documentation updates

  * some fixes to the numerical Jacobian involving X vs. Y

  * a new subCh network for He burning was added.

  * implemented the new c12(a,g)o16 nuclear reaction rate and its
    corresponding inverse from the work of Deboer et al. 2017 
    (https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.89.035007). 
    To use the new rate, user must set `use_c12ag_deboer17` to `true`. 
    This rate is only useable in the `aprox13`, `aprox19`, `aprox21`,
    and `iso7` reaction rate networks. Closes issue #44.

  * a routine util/cj_detonation was added to compute the
     Chapman-Jouguet detonation velocity for any of the networks

  * the burn retry strategy now sticks with the current integrator and 
    uses looser tolerances before switching to a different integrator.

# 18.04

   * pynucastro (https://github.com/pynucastro/pynucastro) can now
     generate reaction networks compatible with StarKiller.  See the
     subch network.

# 17.11

   * a new option to boost the reaction rates has been added
     to the integrators (PR #64)

   * we now disable some composition derivatives in the EOS
     by default, for performance and memory reasons.  They can
     be reenabled by defining the preprocessor variable 
     EXTRA_THERMO (PR #59)

# 17.10

  * the compositional derivatives are no longer available
    by default from the EOS.  To get these, set the preprocessor
    variable EXTRA_THERMO.  This change was done for performance
    reasons.

  * the aprox19 and aprox21 networks no longer use a numerical
    Jacobian by default, as this was found to result in some
    bad numerical issues in VODE (PR #49)

  * the maximum temperature for reactions, MAX_TEMP, is now
    an adjustable input parameter rather than being hardcoded
    at 1.d11.

  * the Helmholtz EOS table is now read by the IO processor and
    broadcast to other processors (PR #53)


  * the VODE integrator now does some additional checks on the
    state to ensure consistency (PR #47)

# 17.09

  * a new rety mechanism was implemented that allows a different
    integrator to be used if the primary integrator fails

  * the electron Ni56 electron capture rates and energy losses
    were updated from Mazurek (1973) to LMP (2000).  Thanks to
    Carl Fields for this contribution.  Pull request #40


# 17.08

  * fix to aprox21 from Aron Michel (HITS) that fills in missing
    reactions

  * updated the helmholtz EOS to use the latest table from Frank
    Timmes (in particular, this is now denser, with 2x points in T and
    rho dimensions).  If you copied the old table, you need to make sure
    you are using the new table now.

  * add stellar conductivities from Frank Timmes


# 17.06

  * a new Fortran 90 port of VODE has been added

  * the unit tests now require AMReX instead of BoxLib to build


# 17.01

  * we've removed the option to integrate molar fractions and instead
    the ODE system always operates on mass fractions (the networks
    return derivatives of molar fractions and they are automatically
    converted).


# 16.12

  * a new unit test, test_sdc, was created to test the SDC interface
    to the networks

  * we now rely on the network module to provide aion_inv (1/aion)

  * the VODE integrator now supports SDC integration


# 16.09

  * num_rate_groups is now a property of the individual networks

  * a new integration method, Rosenbrock, was added to the BS
    option (set ode_method)

  * the number of RHS and Jac evaluations is now passed out
    of the burner through the burn_t type for diagnostic and
    load-balancing use

  * support for spectral deferred correction coupling of the
    burner and hydro was added to the BS integrator


# 16.08

  * Microphysics/eos/ has been renamed Microphysics/EOS/ to better
    conform to the conventions used in Castro and Maestro

  * the User's Guide has been extensively updated

  * OpenMP and OpenACC directives have been added to the unit tests

  * the BS integrator's type, bs_t, has now contains a burn_t
    internally, simplifying the conversion from bs_t to call the
    actual_rhs/jac

  * the rates() component of burn_t was removed.  We no longer
    rely on rate caching

  * we now store the simulation time at the start of the burn as
    t0 in the rpar storage to use as an offset to the integration
    time

  * the species derivatives (dh/dX and de/dX) and enthalpy were
    removed from the burn_t

  * a new option to integrate of X instead of Y was added
    (integrate_molar_fraction = F)

  * integration of networks with nspec_evolve < nspec were fixed
    to now apply the algrebic constraint relating mass fractions
    through a new update_unevolved_species() function

  * the electron capture rate on Ni56 used by aprox19 and aprox21 was
    fixed

  * the BS integrator can now use the initial timestep estimation
    algorithm that VODE uses 9use_timestep_estimator = T)

  * a centered difference numerical Jacobian option was added


# 16.07

  * we now use MICROPHYSICS_HOME instead of MICROPHYSICS_DIR as the
    environment variable to point to the Microphysics/ directory.

  * there are now two standalone unit tests, test_react and test_eos
    that don't need Maestro or Castro to compile.

  * new burn mode that limits numerically unstable burning.

  * UsersGuide/ was renamed to Docs/ to be consistent with the other
    BoxLib codes

  * the energy equation now uses an offset to help with the BS ODE
    integration convergence

  * the runtime parameter small_x now is owned by the network
