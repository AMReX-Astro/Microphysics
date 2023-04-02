# 23.04

  * added preliminary CMake support (#1151, #1164, #1166)

  * added code of conduct (#1152)

  * clang-tidy code clean-ups(#1141, #1153, #1156)

  * burn_t now stores whether we entered NSE (#1144, #1147)

  * burn_t now store chemical potentials for NSE (#1149)

  * some NSE solver updates to make it easier to enter NSE (#1138, #1139)

  * a new CNO + rp network for XRBs (#1145)

# 23.03

  * updated all of the pynucastro networks to the latest
    pynucastro version (#1134, #1136)

  * added tricubic interpolation for the NSE table (#1114)

  * fixed an issue with rate tabulation in the aprox nets (#1123,
    #1124)

  * fixed some bugs in the NSE solver and made the hybrid Powell
    solver more robust (#1122)

# 23.02

  * `T_fixed` now works with NSE (#1098, #1111)

  * `USE_NSE` was changed to `USE_NSE_TABLE` (#1108)

# 23.01

  * a new test, burn_cell_primordial_chem/, works with the primordial
    chemistry (#1064)

  * burn_cell and burn_cell_sdc now work with aux data with NSE
    (#1084, #1094)

  * a screening implementation from Chabrier & Potekhin (1998) was
    added (#1085)

  * test_react can now output the burn state that took the longest to
    evaluate (#967)

  * an initial implementation of adaptive nuclear statistic equilibrium
    was added (#983)

# 22.12

  * A first order backward Euler integrator was added that works with
    both Strang and simplified-SDC integration (#504, #1041, #1042, #1045)

  * The jacobian = 3 option for simplified SDC was removed (#1069)

  * An option to not subtract off the initial energy after the burn
    was added as well as one to evolve number densities (#999, #1051)

  * The python bindings have been removed (#1036)

  * An issue with reading the helmholtz table on GPUs was fixed (#1020)

# 22.11

  * use of the auxillary state to define composition is now enabled
    via USE_AUX_THERMO and the preprocessor variable AUX_THERMO
    (#1003)

# 22.10

  * A null screening routine was added to disable screening for any network at
    compile time. (#992)

  * An option to disable the clipping of species in the VODE integration
    was added (`integrator.do_species_clip`) (#989)

  * A unit test for C++ partition functions was added (#980)

  * An EOS for primordial chemistry was added (#981)

  * A C++ version of the Powell's hybrid non linear system solver was
    added (#976)

  * screening in the approximate rates in pynucastro nets was fixed
    (#978)

# 22.09

  * An NSE solver was added (#963)

  * A new network, subch_simple, was added that further simplifies
    He/C burning (#964)

# 22.08

  * The subch network was renamed `subch_approx` and the subch2
    network was renamed `subch_full` (#947)

# 22.07

  * Two new screening formulations have been added for reaction rates,
    based on Chugunov, DeWitt, and Yakovlev 2007 and Chugunov and
    DeWitt 2009.  These can be used with any network by setting
    SCREEN_METHOD at compile time.(#887)

# 22.06

  * The subch2 network now has runtime parameters allowing for
    some key rates to be disabled (#921).

# 22.05

  * The subch2 network was improved by adding some missing C+C, C+O,
    and O+O rates. (#915)

# 22.04

  * aprox networks now use a templated C++ righthand side formulation
    that builds the ODE system at compile time. (#802)

  * pynucastro networks were regenerated to take advantage of recent
    optimizations (#901)

# 22.02

  * The Microphysics repo was moved to the AMReX-Astro github
    organization: https://github.com/amrex-astro/Microphysics

    You can update your git remote via:

    git remote set-url origin git@github.com:amrex-astro/Microphysics.git

  * Fortran support has been removed from the runtime parameter
    scripts (#869)

# 22.01

  * we added back in support for the "Nonaka plot".  This outputs the
    state in the RHS routine for a single zone during the reaction
    network integration. (#830)

  * we removed the xrb_simple network.  This was never used in any
    science calculations (#827)

  * the simplified-SDC step rejection logic in VODE was improved (#818)

# 21.12

  * all of the pynucastro networks were regenerated with the latest
    pynucastro and converted to C++.  Performance was also improved
    (#809)

  * a bug was fixed in the VODE step rejection logic (#815)

  * we added USE_MICROPHYSICS_DEBUG that defines MICROPHYSICS_DEBUG to
    turn on more verbosity to help with debugging (#817)

# 21.11

  * burn_cell was not correctly doing substepping in some cases.
    This has been fixed (#784)

  * With Intel compilers, logical runtime parameters in Fortran
    were not being correctly cast to int (#789)

  * Simplified-SDC now works with Fortran nets (#786)

# 21.09

  * Added a new nova network (nova2) with pp and (hot-)CNO and some
    breakout reactions (#751)

  * Some fixes to the NSE bailout in aprox19 (#739, #753, #755) and
    the relaxation check on the NSE critera (#754)

  * Added a new unit test for single-zone SDC (burn_cell_sdc) (#744)

# 21.08

  * test_react can now be run with a uniform composition to test GPU
    performance (#734)

  * the numerical Jacobian now uses a more robust one-sided difference
    algorithm (#660, #728)

  * for simplified SDC, we now only integrate (rho X, rho e), and no
    longer integrate (rho E) (#710, #712, #717)

  * for the NSE bailout, we can now relax the conditions needed to enter
    NSE after a failed burn (#702)

# 21.07

   * The C++ networks now implement abort_on_failure functionality
     (#697)

# 21.06

   * The ability to use a system BLAS library was removed (#675)

   * An equation of state for hypervelocity impacts was added (#645)

# 21.05

   * For aprox19 + NSE, we now "bail out" of the integration
     immediately if the state enters NSE, and then do the rest of the
     update through the NSE table. (#658)

   * The old gamma_law EOS was removed and gamma_law_general was
     renamed gamma_law.  The old gamma_law EOS have a very reduced
     subset of thermodynamic quantities that it computed, for
     efficiency purposes.  This is no longer needed now that we have
     templated the EOSes and have different eos_t data types (#653).

   * Integration for simplified-SDC was interpreting rtol incorrectly.
     This has been fixed (#643)

   * Screening for the 3-alpha reaction in the subch, subch2, and nova
     networks was fixed (#627, #634, #635)

# 21.04

   * We added a new mechanism to recover a failed burn when the state
     tries to enter NSE during the evolution, when using the aprox19 +
     NSE network.  Now it will capture the failure and redo the burn
     if it satisfies the NSE criteria (#628)

   * We updated the VODE logic for rejecting a step to consider mass
     fractions for both simplified-SDC and true SDC burns (#619)

# 21.03

   * We now integrate internal energy (e) directly instead of
     integrating temperature (T) for the thermodynamic evolution. T is
     obtained from e with an EOS call when needed to evaluate the
     rates. (#496)

   * simplified-SDC can be used with the NSE table in aprox19 now
     (#423, #497)

# 21.02

   * Fortran support for the VODE integrator has been removed (#538)

   * Runtime parameters can now be set in the inputs file instead of
     the probin file (and then they are read in by C++ ParmParse).  If
     a parameter is set in both places, then the inputs value is used.
     (#505)

   * Fortran support for simplified-SDC in the VODE integrator has
     been removed. (#492)

# 21.01

   * Microphysics now requires C++17 (gcc >= 7, CUDA >= 11). (#485)

   * The BS integrator was removed.  This was Fortran only, doesn't
     support SDC integration, and not well used. (#488)

# 20.12

   * The default absolute tolerance for species (atol_spec) has been
     increased to 1.e-8 (from 1.e-12). (#422)

   * An interface has been added for C++ integrators to call the RHS
     from a network that only has a Fortran implementation. This allows
     the use of USE_CXX_REACTIONS = TRUE for any network (however, CUDA
     is not currently supported for this case). (#419)

# 20.11

   * The aprox19 + NSE network was ported to C++ (#362)

   * The simplified-SDC code path was ported to C++ (#389)

# 20.10

   * An option to use NSE instead of integrating the reaction
     network has been added to the aprox19 network. (#332)

   * The BS integrator no longer supports simplified-SDC (#393)

   * The triple_alpha_plus_cago network switch to using binding
     energies in MeV, consistent with the aprox nets (#354)

# 20.09

   * Unit tests now write a job_info file (#383)

   * A new single-zone EOS test routine was created as unit_test/eos_cell
     (#382)

   * The gamma_law eos (not gamma_law_general) now fills the sound
     speed, entropy, and derivatives for more inputs (#374)

   * The rprox network now has screening (#377)

   * The NETWORK_PROPERTIES file was split to put the number of
     auxiliary species into its own file, NAUX_NETWORK.  This allows
     us to put if-logic into the file to choose the number of
     auxiliary quantities based on make setting (like USE_NSE).
     (#370)

# 20.08

   * Several of the unit tests had separate C++ and Fortran
     implementions.  These have been unified (#343, #344, #345)

   * The VBDF integrator was removed (#348)

   * VODE can now reject an internal timestep that has any abundance
     change by more than a factor of 2, or an abundance < 0 or > 1, as
     well as timesteps where the temperature ends up negative. (#350)

# 20.07

   * The "master" branch has been renamed "main" (#333)

   * NETWORK_PROPERTIES now includes the number of Aux quantities (#330)

# 20.06

   * For integration with simplified SDC, we now interpret atol_spec
     as an absolute tolerance on X alone instead of (rho X) (#311)

   * burn_cell can now use the C++ burner if compiled with
     USE_CXX_REACTIONS=TRUE and run with do_cxx = 1. (#313)

   * The original burn_cell (which used the F90 BoxLib build system)
     is removed and replaced with burn_cell_C (which uses the newer
     build system). (#316)

   * The analytic Jacobian with simplified SDC now is written in terms
     of the conserved fluid state and works for a wide range of
     problems (#228)

# 20.05

   * We now have an option for using sparse storage for aprox13 in C++
     (#307)

   * iso7 and aprox13 are now available as a C++ network (#303, 305)

   * species names are available as an enum in network_properties.H (#304)

   * The screening on O16+O16 in iso7 was fixed (#302)

   * The VODE integrator is now available in C++ (#299)

# 20.04

   * The wion network property was removed (#294)

   * There are new unit tests for the screening and aprox rates
     modules (both C++ and Fortran interfaces).

   * The screening routines were ported to C++ (#290) and the screenz
     routine was removed in favor of screen5 (#293)

   * a new method, is_input_valid, was added to all EOSes (both C++
     and Fortran interfaces) that can be used to query whether an EOS
     supports a particular input mode (e.g. eos_input_rp).  (#291)

   * The aprox rates used with iso7, aprox13, aprox19, and aprox21
     have been converted to C++ (#288)

   * We've rewritten the VODE integrator to remove all "go to"
     statements (#275, 276, 278, 280, 281, 282, 283, 284, 285, 286,
     287)

   * We removed the ability to have nspec_evolve < nspec.  This
     feature was not widely used and greatly complicated the code
     paths.(#279)


# 20.03

   * The nuclei information for both Fortran and C++ is now
     automatically generated from a network inputs file at compile
     time.  As part of this change, 1/A is precomputed and stored as a
     constant (#253, 258).

   * The license for StarKiller Microphyscs was made explicit and
     a license.txt file was added (#267)

   * A framework for pure C++ EOSes has been created and a pure C++
     unit test, test_eos_C, is available to test these.  (#246) The
     following EOSes have been ported to C++: ztwd (#268), multigamma
     (#265), polytrope (#264), gamma_law (#263), helmholtz (#262),
     gamma_law_general (#246), rad_power_law (#269), breakout (#270)

   * The GPackage.mak files that were a remnant of the old
     BoxLib F90 build system have been removed.  They were
     not maintained.  (#212).

   * All of the interface files have been collected together
     in the interfaces/ dir.  (#240)

   * The network C++ headers have been renamed network_properties.H
     and the nuclei information (aion and zion) have been
     added. (#244)

# 20.02

   * Added a C++ header file, actual_network.H, that defines the
     network size.  This is the start of making making the
     microphysics routines available in C++.

   * regenerated the pynucastro networks with the latest weak rate
     formulations from pynucastro.

# 20.01

   * The burn_t type no longer includes ydot or jac -- this allows
     us to optimize the memory access on GPUs (#220)

   * The radiation pressure contribution to the Helmholtz EOS has
     had a dampener applied to it that makes it approximately zero
     for low densities where the radiation pressure would lead to
     unphysical situations like a superluminal sound speed. (#235)

   * The original VODE integrator was removed and the Fortran 90
     version VODE90 was renamed to VODE. (#221)

   * The test_react unit tests no longer require a composition inputs
     file (xin*).  They now create the composition profile at runtime.
     (#211)

# 19.12

   * Simplified SDC integration now uses the same retry strategy
     as the default (non-SDC) integration. (#215)

   * VODE90 can now participate in the retry strategy that was
     previously available to the VODE integrator, where it can
     switch to the BS integrator if loosening the tolerances
     does not allow the burn to complete. (#201)

   * The parameter ode_max_steps was made consistent in VODE and
     VODE90; in some places it was being ignored. (#214)

   * The helmholtz EOS was restructured, splitting the different
     components into different functions and optimizing the memory
     accesses. (#200)

   * The derivatives with respect to mass fraction (dpdX, dedX, dhdX)
     were removed from eos_t and are now available through a new type,
     eos_xderivs_t and the composition_derivatives() routine.  (#207)

   * A bug in the screening of the C12+C12 and O16+O16 rates in iso7
     was fixed.

   * The test_eos unit test now outputs all of the variables in the
     eos_t type.

# 19.11

   * VODE90 now works with the simplified SDC time step algorithms,
     and the preprocessor option for this SDC was changed to
     SIMPLIFIED_SDC (#194)

   * rprox now works on GPUs

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
