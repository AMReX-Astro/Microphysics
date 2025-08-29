# Changelog

## 25.08

  * work on the JOSS paper (#1838, #1839)

  * remove old inputs files from pynucastro networks (#1836)

  * updates to the NSE solver (#1829, #1833)

  * fix string_view in some physics (#1832)

  * new CI for validating JSON (#1831)

## 25.07

  * module string names are now `constexpr std::string_view` (#1825)

  * doc updates (#1796, #1809, #1810, #1818, #1824, #1827, #1830)

  * JOSS paper (#1658)

  * the ASE network has been updated to better support NSE (#1813)

  * `hybrj` now can take the Jacobian and constrain equations as input
    (#1817)

  * remove amrex namespace (#1820)

  * for He-burning nets, don't consider Suzuki rates in pynucastro
    (#1815)

  * `nse_sdc_burn` has been updated for `NSE_NET` to give better NSE
    agreement with the networks (#1812)

  * `sneut5` now used autodiff (#1799, #1808)

  * recombination neutrinos are now disabled by default in `sneut5` (#1793, #1794)

  * the `nse_compatibility` script has been updated (#1811)

  * autodiff improvements and optimizations (#1803)

  * enable more `clang-tidy` checks (#1807)

  * clean up the `nova` network script (#1804)

  * regenerate the non-He burning pynucastro nets (#1797)

  * update `sn160` to include tabular rates (#1805)

  * update the constants to use the same CODATA values as pynucastro
    (#1802)

  * allow constant temperature in the SDC burn (#1801)

## 25.06

  * update `README.md` (#1785)

  * doc updates (#1787, #1795)

  * `test_eos` now uses the composition runtime parameters (#1792)

  * code cleaning (#1790, #1791)

  * neutrino cooling can now has an option to remove the recombination
    contribution (#1789)

## 25.05

  * some clang-tidy cleaning (#1784)

  * with HIP we were disabling inlining due to ROCm issues.  This is
    now an option, with the default enabling inlining, since it works
    as expected for ROCm >= 6.3.1 (#1780)

  * clean up the pynucastro scripts that derived reverse rates (#1781)

## 25.04

  * the "he-burn" networks now will look for weak rates for
    all nuclei, not just the Fe-group (#1763)

  * clean up the ECSN network generation script (#1776)

  * Allow for single-step backward Euler (#1773)

## 25.03

  * the `nova2` net was renamed `nova`, and the old `nova`
    net was removed (#1746, #1768)

  * small improvements to the he-burn scripts (#1765)

  * documentation updates (#1741, #1744, #1747, #1751, #1756, #1759,
    #1760, #1766, #1767)

  * restrict Sphinx < 8.2.0 due to nbsphinx issues (#1762)

  * centralize some common unit test runtime parameters (#1749) and
    explicitly set `small_dens` and `small_temp` for some tests
    (#1757)

  * remove old, unneeded macros from the build system (#1717, #1753,
    #1754, #1755)

  * modernize parallel loops in some unit tests (#1752)

  * update the HIP CI action to a later runner

  * namespace and clang-tidy fixes (#1745)

## 25.02

  * documentation updates (#1678, #1690, #1692, #1695, #1700, #1702,
    #1703, #1707, #1709, #1712, #1713, #1714, #1715, #1716, #1723,
    #1726, #1732, #1733) including renaming the docs directory to
    `Docs` (#1689) and the addition of a CITATION.md (#1731)

  * codespell fixes (#1724)

  * rename `test_nse_net` -> `nse_net_cell`; `test_part_func` ->
    `part_func_cell` (#1729)

  * remove old testing scripts (#1727)

  * update `test_react` `README.md` (#1722)

  * implement Debye-Huckel screening and allow it to be used as a test
    for whether screening is needed (#1688)

  * remove the `use_raw_inputs` option from the EOS (#1721)

  * remove some old Fortran references (#1718, #1719)

  * switch from `std::clamp` to `amrex::Clamp` due to HIP compiler
    issues (#1711)

  * reorganize the He nets -- they are now all under
    `networks/he-burn` and share a common python setup (#1687, #1710)
    also update these nets with pynucastro (#1685).  This also
    add a new network with 31.

    Renamed nets are:
    * `subch_base` -> `he-burn/he-burn-18a`
    * `subch_simple` -> `he-burn/he-burn-22a`
    * `He-C-Fe-group` -> `he-burn/he-burn-36a`
    * `CNO_He_burn` -> `he-burn/cno-he-burn-33a`

  * fix the Chapman-Jouguet detonation utility (#1699)

  * update the mailmap (#1706)

  * switch `std::pow(x, 1./3.)` to `std::cbrt(x)` (#1705)

  * remove `do_acc` option (#1708)

  * make the breakout EOS check if 1/mu is defined (#1694)

  * remove old Doxygen (#1697) and skynet (#1698) scripts from `util/`

## 25.01

  * update HIP/CUDA dependences to include sparse libraries (#1686)

  * rename `Opacity_dir` -> `OPACITY_DIR` (#1679)

  * update the integration and NSE docs (#1682)

## 24.12

  * documentation improvements (#1661, #1667, #1670)

  * optimize tabular NSE EOS calls (#1668)

  * CI fixes (#1666, #1671, #1675) and new partition function CI
    (#1673)

  * `burn_cell` can now initialize all mass fractions to be equal
    (#1665)

## 24.10

  * metal chemistry updates (#1648) with ices (#1650) and cosmic rays (#1651)

  * added dust to primordial chemistry (#1649)

  * doc updates (#1652)

## 24.09

  * Improvements to the primordial chemistry network and the addition
    of a new version that includes metals and dust (#1642, #1644)

  * code clean-ups (#1645)

  * documentation improvements (#1637)

  * outputting the `burn_t` now prints the mass fractions / number densities
    in scientific notation (#1643)

  * improvements to the looping and zeroing of the Jacobian in the
    integrators (#1636, #1640)

## 24.08

  * autodiff is now used with the templated reaction networks (#1614)
    and some autodiff clean-ups and derivative fixes (#1604, #1612,
    #1613, #1616, #1619, #1633)

  * we can now output warnings from GPUs if you compile with
    `USE_GPU_PRINTF=TRUE` (#1629, #1635)

  * documentation improvements (#1570, #1628)

  * a new jacobian unit (`jac_cell`) test was added that compares the
    numerical and analytic Jacobians (#1618)

  * support for Strang + NSE has been removed.  NSE only works with
    SDC now (#1549, #1621)

  * the network `CNO_He_burn` was added for explosive H/He burning
    (#1622)

  * code clean-ups (#1582, #1602, #1609, #1623, #1624, #1625, #1626,
    #1627, #1631, #1639)

  * `test_nse_net` now also tests the NSE EOS interface (#1621)

  * the self-consistent NSE + SDC update has been synced with the
    tabular NSE implementation (#1569, #1607, #1617)

  * `test_jac` was not correctly evaluating the numerical Jacobian
    (#1615)

  * the `fast_atan` function is now more accurate (#1611)

  * `test_ase` was renamed `test_nse_net` and the old `test_nse` was
    removed (#1610)

  * the old `test_screening` unit test was removed (#1608)

  * the RKC integrator now supports NSE bailout (#1544)

  * a second temperature check for tabular NSE was added -- above this
    temperature we don't consider composition (#1547)

  * a SDC+NSE unit test was added (#1548)

  * a fast log and fast pow approximation was added (#1591)

  * the `primordial_chem network` now uses the fast math routines (#1605)

  * fix potential Inf in constexpr linear algebra (#1603)

## 24.07

   * added an autodiff library and converted all of the screening
     functions to use autodiff for the thermodynamic derivatives
     (#1581, #1588, #1593, #1596, #1597, #1600)

   * some testing infrastructure fixes (#1598, #1589)

   * documentation improvements (#1594)

   * added approximate math functions for exp and atan (#1583, #1586)

   * fix return code for PrimordialChem unit test (#1590)

   * NSE optimizations (including `chabrier1998` screening) (#1585)

   * remove "using namespace amrex" from most headers (#1584)

   * NSE table can work with other network now (#1576, #1577, #1580)

   * the `subch_full` and `subch_approx` networks were removed -- these
     are replaced by `subch_simple` and `subch_base` (#1578)

   * retry tolerances now default to use the same values as the first
     attempt, unless they are explicitly set in an inputs file (#1573)

## 24.06

   * added the ability to access the runtime parameters via a struct.
     This will eventually be used to remove the dependency on globals
     (#1433, #1575)

   * simplified the integrators by extracting common infrastructure
     into setup and cleanup functions (#1515, #1546)

   * lots of documentation improvements including sectioning (#1559)
     integrators (#1560, #1561, #1562, #1563, #1566, #1567, #1568),
     runtime parameters (#1557), and link checking (#1552)

   * CUDA no longer requires VODE + generalize some `AMREX_USE_CUDA`
     to `AMREX_USE_GPU` (#1564)

   * self-consistent NSE now accepted the temperature threshold as a
     runtime parameter (#1558)

   * general code cleanups (#1537, #1551, #1553, #1554)

   * unit tests no longer append `.cxx` to output (#1309)

   * added an `eos_rh_t` EOS type (#1539)

## 24.05

   * Runtime parameters can now be type `bool` (#1536)

   * more clang-tidy and compiler warning cleaning (#1527 #1530,
     #1532, #1533)

   * Remove recursion in quicksort to avoid CUDA stack limits (#1531)

   * Update the pynucastro networks to cache derived rate partition
     functions (#1529)

## 24.04

   * A new `test_screening_templated` unit test was added -- this
     works with any of the templated networks. (#1525)

   * A lot of small code clean-ups from clang-tidy (#1516, #1518, #1519, #1520, #1522)

   * The NSE solver was optimized (#1503, #1504, #1506, #1507, #1508)

   * The integrator code was synced up between implementations, fixing
     a bug in the RKC retry tolerances (#1513)

   * A `reinterpret_cast` in `rhs.H` was removed (#1435)

## 24.03

   * pivoting in the linear algebra routines can now be disabled
     (#1454)

   * the scaling of the energy derivatives in the Jacobian when
     running with `integrator.scale_system=1` has been fixed (#1479)

   * added a new linear algebra unit test (#1493)

   * when building with HIP we disable forced inlining (#1490)

   * improved the energy update with NSE and remove unused terms
     (#1483, #1484, #1485)

   * remove `using namespace amrex` from most headers (#1465, #1474)

   * updated the pynucastro networks to pynucastro 2.2.0 (#1470)

   * fixed an error code check in the VODE integrator (#1472)

   * added a zone-by-zone retry capability to the burner (#969)

## 24.02

   * Lots of general code cleaning from coverity and clang-tidy
     (#1450, #1452, #1453, #1460, #1459, #1457, #1458)

   * Fixed a bug in the VODE pivoting when a cached Jacobian is used
     (#1456)

   * Added reverse rates to `CNO_extras` (#1445)

   * Sync networks up with pynucastro to get `constexpr` mion/bion
     (#1437)

   * NSE+SDC improvements (#1431)

   * Start of moving the runtime parameters from globals to structs
     (#1422)

## 24.01

   * The quantum corrections for the Chabrier screening are
     now optional (#1428)

   * We've replaced `std::pow()` with `amrex::Math::powi` for integer
     powers for better GPU performance (#1432)

   * `in_nse` now works with an `eos_t` for tabular NSE (#1424)

   * There are a new set of interfaces for inverting the EOS when
     we are in NSE (with the tabular NSE) that consistently find
     T and the composition (#1405, #1430)

   * The NSE table now uses finer spacing (#1427)

   * The SDC+NSE update for tabular NSE is now based on a 2nd-order
     Runge-Kutta method (#1415)

   * An additional check on molar fractions was added to
     self-consistent NSE to determine if we are in NSE (#1417)

   * The NSE table interface was changed (#1404, #1418)

   * A script that checks if a network is compatible with
     self-consistent NSE was added (#1416)

   * constant T evolution was fixed (#1408)

   * An unsafe `reinterpret_cast` was removed from linear algebra (#1412)

   * Methods for computing T derivatives of an NSE table quantity were
     added (#1407)

## 23.12

  * The SDC+NSE update now includes plasma neutrino losses (#1357,
    #1400)

  * The default tabular NSE interpolation is now cubic (#1399)

  * Self-consistent NSE now requires `chabrier1998` or `null` screening
    (#1398)

  * A new network, `subch_base`, was added that further simplifies
    `subch_simple` (#1393)

  * A slightly larger network for Urca was added (#1365)

  * A new NSE table was added.  This is generated via pynucastro and
    there is a python script that can be used to regenerate it (#1350)

  * A bug was fixed in the neutrino cooling that was introduced in an
    optimization last release (#1380)

## 23.11

  * The `sneut5` neutrino cooling term was cleaned up (#1371, #1372,
    #1373, #1374, #1375, #1377, #1378, #1379)

  * The number of predictor-corrector iterations for the SDC+NSE algorithm
    is now a runtime parameter (#1370)

  * The Urca network now includes a more accurate rate for neutron decay
    and electon-capture onto a proton. (#1359)

  * The `He-C-Fe-group` network now includes the positron parts of the
    weak reaction rates (#1360)

  * A check was added to ensure that the `helm_table.dat` is valid on
    reading (#1355)

## 23.10

  * The simplified-SDC and true-SDC code paths for integration
    have been merged (#1338, #1340, #1341).

  * All pynucastro networks have been updated with the latest
    version of pynucastro (2.1.0) (#1342)

  * The neutrino cooling terms now use templating on derivatives
    (#1329)

  * `NUM_EXTRA_SPECIES` was removed (#1321)

## 23.09

  * The file `NETWORK_PROPERTIES` has been removed from each network,
    as the legacy screening method is no longer used. (#1310)

  * The `rprox` network was updated and the Jacobian was fixed (#1300)

  * The `primordial_chem` EOS now can take density and pressure as
    inputs (#1302)

## 23.07

  * The preprocessor variable `EXTRA_THERMO` has been removed.
    Use cases that depend on dpdA/dpdZ or dedA/dedZ should use
    `eos_extra_t`, which is a container that holds all of the
    entities in `eos_t` as well as these derivatives wrt A and Z. (#1229)

  * added the ability to scale the energy we integrate by
    the initial energy in the ODE integration (#1224)

  * added an implementation of the Gershgorin circle theorem
    for estimating the spectral radius of our ODE system (#1222)

  * removed `SDC_EVOLVE_ENTHALPY` -- this was not being used (#1204)

## 23.06

  * Added a new Runge-Kutta-Chebyshev integrator (#1191)

  * Lots of clean-up to the primordial chem network (#1180, #1181
    #1198)

  * Namespaces for the runtime parameters are now required in C++ (#1056)

  * The SDC+NSE update for tabular NSE was fixed -- we were previously
    computing the energy release incorrectly (#1092)

## 23.05

  * The `abort_on_failure` runtime parameter has been removed (#1174)

## 23.04

  * added preliminary CMake support (#1151, #1164, #1166)

  * added code of conduct (#1152)

  * clang-tidy code clean-ups(#1141, #1153, #1156)

  * `burn_t` now stores whether we entered NSE (#1144, #1147)

  * `burn_t` now store chemical potentials for NSE (#1149)

  * some NSE solver updates to make it easier to enter NSE (#1138, #1139)

  * a new CNO + rp network for XRBs (#1145)

## 23.03

  * updated all of the pynucastro networks to the latest
    pynucastro version (#1134, #1136)

  * added tricubic interpolation for the NSE table (#1114)

  * fixed an issue with rate tabulation in the aprox nets (#1123,
    #1124)

  * fixed some bugs in the NSE solver and made the hybrid Powell
    solver more robust (#1122)

## 23.02

  * `T_fixed` now works with NSE (#1098, #1111)

  * `USE_NSE` was changed to `USE_NSE_TABLE` (#1108)

## 23.01

  * a new test, `burn_cell_primordial_chem`, works with the primordial
    chemistry (#1064)

  * `burn_cell` and `burn_cell_sdc` now work with `aux` data with NSE
    (#1084, #1094)

  * a screening implementation from Chabrier & Potekhin (1998) was
    added (#1085)

  * `test_react` can now output the burn state that took the longest to
    evaluate (#967)

  * an initial implementation of adaptive nuclear statistic equilibrium
    was added (#983)

## 22.12

  * A first order backward Euler integrator was added that works with
    both Strang and simplified-SDC integration (#504, #1041, #1042, #1045)

  * The jacobian = 3 option for simplified SDC was removed (#1069)

  * An option to not subtract off the initial energy after the burn
    was added as well as one to evolve number densities (#999, #1051)

  * The python bindings have been removed (#1036)

  * An issue with reading the helmholtz table on GPUs was fixed (#1020)

## 22.11

  * use of the auxiliary state to define composition is now enabled
    via `USE_AUX_THERMO` and the preprocessor variable `AUX_THERMO`
    (#1003)

## 22.10

  * A `null` screening routine was added to disable screening for any network at
    compile time. (#992)

  * An option to disable the clipping of species in the VODE integration
    was added (`integrator.do_species_clip`) (#989)

  * A unit test for C++ partition functions was added (#980)

  * An EOS for primordial chemistry was added (#981)

  * A C++ version of the Powell's hybrid non linear system solver was
    added (#976)

  * screening in the approximate rates in pynucastro nets was fixed
    (#978)

## 22.09

  * An NSE solver was added (#963)

  * A new network, `subch_simple`, was added that further simplifies
    He/C burning (#964)

## 22.08

  * The `subch` network was renamed `subch_approx` and the `subch2`
    network was renamed `subch_full` (#947)

## 22.07

  * Two new screening formulations have been added for reaction rates,
    based on Chugunov, DeWitt, and Yakovlev 2007 and Chugunov and
    DeWitt 2009.  These can be used with any network by setting
    `SCREEN_METHOD` at compile time. (#887)

## 22.06

  * The `subch2` network now has runtime parameters allowing for
    some key rates to be disabled (#921).

## 22.05

  * The `subch2` network was improved by adding some missing C+C, C+O,
    and O+O rates. (#915)

## 22.04

  * aprox networks now use a templated C++ righthand side formulation
    that builds the ODE system at compile time. (#802)

  * pynucastro networks were regenerated to take advantage of recent
    optimizations (#901)

## 22.02

  * The Microphysics repo was moved to the AMReX-Astro github
    organization: https://github.com/amrex-astro/Microphysics

    You can update your git remote via:

    git remote set-url origin git@github.com:amrex-astro/Microphysics.git

  * Fortran support has been removed from the runtime parameter
    scripts (#869)

## 22.01

  * we added back in support for the "Nonaka plot".  This outputs the
    state in the RHS routine for a single zone during the reaction
    network integration. (#830)

  * we removed the `xrb_simple` network.  This was never used in any
    science calculations (#827)

  * the simplified-SDC step rejection logic in VODE was improved (#818)

## 21.12

  * all of the pynucastro networks were regenerated with the latest
    pynucastro and converted to C++.  Performance was also improved
    (#809)

  * a bug was fixed in the VODE step rejection logic (#815)

  * we added `USE_MICROPHYSICS_DEBUG` that defines `MICROPHYSICS_DEBUG` to
    turn on more verbosity to help with debugging (#817)

## 21.11

  * `burn_cell` was not correctly doing substepping in some cases.
    This has been fixed (#784)

  * With Intel compilers, `logical` runtime parameters in Fortran
    were not being correctly cast to `int` (#789)

  * Simplified-SDC now works with Fortran nets (#786)

## 21.09

  * Added a new nova network (`nova2`) with pp and (hot-)CNO and some
    breakout reactions (#751)

  * Some fixes to the NSE bailout in `aprox19` (#739, #753, #755) and
    the relaxation check on the NSE criteria (#754)

  * Added a new unit test for single-zone SDC (`burn_cell_sdc`) (#744)

## 21.08

  * `test_react` can now be run with a uniform composition to test GPU
    performance (#734)

  * the numerical Jacobian now uses a more robust one-sided difference
    algorithm (#660, #728)

  * for simplified SDC, we now only integrate (rho X, rho e), and no
    longer integrate (rho E) (#710, #712, #717)

  * for the NSE bailout, we can now relax the conditions needed to enter
    NSE after a failed burn (#702)

## 21.07

   * The C++ networks now implement `abort_on_failure` functionality
     (#697)

## 21.06

   * The ability to use a system BLAS library was removed (#675)

   * An equation of state for hypervelocity impacts was added (#645)

## 21.05

   * For aprox19 + NSE, we now "bail out" of the integration
     immediately if the state enters NSE, and then do the rest of the
     update through the NSE table. (#658)

   * The old `gamma_law` EOS was removed and `gamma_law_general` was
     renamed `gamma_law`.  The old `gamma_law` EOS have a very reduced
     subset of thermodynamic quantities that it computed, for
     efficiency purposes.  This is no longer needed now that we have
     templated the EOSes and have different `eos_t` data types (#653).

   * Integration for simplified-SDC was interpreting rtol incorrectly.
     This has been fixed (#643)

   * Screening for the 3-alpha reaction in the `subch`, `subch2`, and `nova`
     networks was fixed (#627, #634, #635)

## 21.04

   * We added a new mechanism to recover a failed burn when the state
     tries to enter NSE during the evolution, when using the `aprox19` +
     NSE network.  Now it will capture the failure and redo the burn
     if it satisfies the NSE criteria (#628)

   * We updated the VODE logic for rejecting a step to consider mass
     fractions for both simplified-SDC and true SDC burns (#619)

## 21.03

   * We now integrate internal energy (e) directly instead of
     integrating temperature (T) for the thermodynamic evolution. T is
     obtained from e with an EOS call when needed to evaluate the
     rates. (#496)

   * simplified-SDC can be used with the NSE table in aprox19 now
     (#423, #497)

## 21.02

   * Fortran support for the VODE integrator has been removed (#538)

   * Runtime parameters can now be set in the inputs file instead of
     the probin file (and then they are read in by C++ ParmParse).  If
     a parameter is set in both places, then the inputs value is used.
     (#505)

   * Fortran support for simplified-SDC in the VODE integrator has
     been removed. (#492)

## 21.01

   * Microphysics now requires C++17 (gcc >= 7, CUDA >= 11). (#485)

   * The BS integrator was removed.  This was Fortran only, doesn't
     support SDC integration, and not well used. (#488)

## 20.12

   * The default absolute tolerance for species (`atol_spec`) has been
     increased to `1.e-8` (from `1.e-12`). (#422)

   * An interface has been added for C++ integrators to call the RHS
     from a network that only has a Fortran implementation. This allows
     the use of `USE_CXX_REACTIONS = TRUE` for any network (however, CUDA
     is not currently supported for this case). (#419)

## 20.11

   * The `aprox19` + NSE network was ported to C++ (#362)

   * The simplified-SDC code path was ported to C++ (#389)

## 20.10

   * An option to use NSE instead of integrating the reaction
     network has been added to the `aprox19` network. (#332)

   * The BS integrator no longer supports simplified-SDC (#393)

   * The `triple_alpha_plus_cago` network switch to using binding
     energies in MeV, consistent with the aprox nets (#354)

## 20.09

   * Unit tests now write a `job_info` file (#383)

   * A new single-zone EOS test routine was created as `unit_test/eos_cell`
     (#382)

   * The `gamma_law` eos (not `gamma_law_general`) now fills the sound
     speed, entropy, and derivatives for more inputs (#374)

   * The `rprox` network now has screening (#377)

   * The `NETWORK_PROPERTIES` file was split to put the number of
     auxiliary species into its own file, `NAUX_NETWORK`.  This allows
     us to put if-logic into the file to choose the number of
     auxiliary quantities based on make setting (like `USE_NSE`).
     (#370)

## 20.08

   * Several of the unit tests had separate C++ and Fortran
     implementations.  These have been unified (#343, #344, #345)

   * The VBDF integrator was removed (#348)

   * VODE can now reject an internal timestep that has any abundance
     change by more than a factor of 2, or an abundance < 0 or > 1, as
     well as timesteps where the temperature ends up negative. (#350)

## 20.07

   * The `master` branch has been renamed `main` (#333)

   * `NETWORK_PROPERTIES` now includes the number of Aux quantities (#330)

## 20.06

   * For integration with simplified SDC, we now interpret `atol_spec`
     as an absolute tolerance on X alone instead of (rho X) (#311)

   * `burn_cell` can now use the C++ burner if compiled with
     `USE_CXX_REACTIONS=TRUE` and run with `do_cxx = 1`. (#313)

   * The original `burn_cell` (which used the F90 BoxLib build system)
     is removed and replaced with `burn_cell_C` (which uses the newer
     build system). (#316)

   * The analytic Jacobian with simplified SDC now is written in terms
     of the conserved fluid state and works for a wide range of
     problems (#228)

## 20.05

   * We now have an option for using sparse storage for `aprox13` in C++
     (#307)

   * `iso7` and `aprox13` are now available as a C++ network (#303, #305)

   * species names are available as an `enum` in `network_properties.H` (#304)

   * The screening on O16+O16 in `iso7` was fixed (#302)

   * The VODE integrator is now available in C++ (#299)

## 20.04

   * The `wion` network property was removed (#294)

   * There are new unit tests for the screening and aprox rates
     modules (both C++ and Fortran interfaces).

   * The screening routines were ported to C++ (#290) and the `screenz`
     routine was removed in favor of `screen5` (#293)

   * a new method, `is_input_valid`, was added to all EOSes (both C++
     and Fortran interfaces) that can be used to query whether an EOS
     supports a particular input mode (e.g. `eos_input_rp`).  (#291)

   * The aprox rates used with `iso7`, `aprox13`, `aprox19`, and `aprox21`
     have been converted to C++ (#288)

   * We've rewritten the VODE integrator to remove all "go to"
     statements (#275, #276, #278, #280, #281, #282, #283, #284, #285,
     #286, #287)

   * We removed the ability to have `nspec_evolve` < `nspec`.  This
     feature was not widely used and greatly complicated the code
     paths. (#279)


## 20.03

   * The nuclei information for both Fortran and C++ is now
     automatically generated from a network inputs file at compile
     time.  As part of this change, 1/A is precomputed and stored as a
     constant (#253, 258).

   * The license for StarKiller Microphyscs was made explicit and
     a `license.txt` file was added (#267)

   * A framework for pure C++ EOSes has been created and a pure C++
     unit test, `test_eos_C`, is available to test these.  (#246) The
     following EOSes have been ported to C++: `ztwd` (#268), `multigamma`
     (#265), `polytrope` (#264), `gamma_law` (#263), `helmholtz` (#262),
     `gamma_law_general` (#246), `rad_power_law` (#269), `breakout` (#270)

   * The `GPackage.mak` files that were a remnant of the old
     BoxLib F90 build system have been removed.  They were
     not maintained.  (#212).

   * All of the interface files have been collected together
     in the `interfaces/` dir.  (#240)

   * The network C++ headers have been renamed `network_properties.H`
     and the nuclei information (aion and zion) have been
     added. (#244)

## 20.02

   * Added a C++ header file, `actual_network.H`, that defines the
     network size.  This is the start of making making the
     microphysics routines available in C++.

   * regenerated the pynucastro networks with the latest weak rate
     formulations from pynucastro.

## 20.01

   * The `burn_t` type no longer includes `ydot` or `jac` -- this allows
     us to optimize the memory access on GPUs (#220)

   * The radiation pressure contribution to the Helmholtz EOS has
     had a dampener applied to it that makes it approximately zero
     for low densities where the radiation pressure would lead to
     unphysical situations like a superluminal sound speed. (#235)

   * The original VODE integrator was removed and the Fortran 90
     version VODE90 was renamed to VODE. (#221)

   * The `test_react` unit tests no longer require a composition inputs
     file (`xin*`).  They now create the composition profile at runtime.
     (#211)

## 19.12

   * Simplified SDC integration now uses the same retry strategy
     as the default (non-SDC) integration. (#215)

   * VODE90 can now participate in the retry strategy that was
     previously available to the VODE integrator, where it can
     switch to the BS integrator if loosening the tolerances
     does not allow the burn to complete. (#201)

   * The parameter `ode_max_steps` was made consistent in VODE and
     VODE90; in some places it was being ignored. (#214)

   * The helmholtz EOS was restructured, splitting the different
     components into different functions and optimizing the memory
     accesses. (#200)

   * The derivatives with respect to mass fraction (dpdX, dedX, dhdX)
     were removed from `eos_t` and are now available through a new type,
     `eos_xderivs_t` and the `composition_derivatives()` routine.  (#207)

   * A bug in the screening of the O16+O16 rate in iso7 was
     fixed. (#204)

   * The `test_eos` unit test now outputs all of the variables in the
     `eos_t` type.

## 19.11

   * VODE90 now works with the simplified SDC time step algorithms,
     and the preprocessor option for this SDC was changed to
     `SIMPLIFIED_SDC` (#194)

   * `rprox` now works on GPUs

## 19.10

   * The `iso7` network was ported to GPUs (#172)

   * VODE90 now better agrees with VODE (#192)

   * When building with multiple integrators, the contents of the `rpar`
     modules could clash.  This has been fixed. (#136)

   * A module for making the "Nonaka plot" tracking the evolution of a
     quantity during the burn was added, and is enabled with
     `USE_NONAKA_PLOT=TRUE` (#190)

## 19.08

   * A new network, `subch2`, was added that combines `aprox13` and the
     `subch` networks. (#184)

## 19.05

   * The `aprox21` network was missing the analytic Jacobian term for
     the derivative of He4 with respect to Ni56. This is fixed. (#175)

   * The numerical Jacobian module, used by the BS and VBDF integrators
     had some wrong scalings.  These have now been fixed (#179, #180)

## 19.01

  * the docs are now automatically build from the sphinx source
    using travis on github.

## 18.12

  * simplify conductivity interface to match the eos interface by
    moving the conductivity into the eos type

## 18.11

  * new python bindings have been added

  * the documentation has been switched to Sphinx and is now hosted
    online.

  * a bug was fixed in the stellarcollapse EOS where the energy offset
    was being applied incorrectly.


## 18.10

  * `test_eos` and `test_react` now both work on GPUs (using the AMReX
    `gpu` branch)

  * the intermediate blending of the weak and strong screening regimes
    was wrong, and has been fixed.  We've also synced some parameters
    up to agree with those in MESA and Kepler.  (#149, #150)

  * `eos_input_is_constant` is now set to `true` for the helmholtz EOS.
    This mean that the EOS inputs will not be modified after the EOS
    call.  This is good for conserving energy in a hydro code, but the
    tradeoff is a small (to root finding tolerance) inconsistency in
    the thermodynamic state. (#154)

## 18.09

  * The Helmholtz parameters `ttol` and `dtol` (controlling the error
    for the Newton iteration when in a mode other than `eos_input_rt`)
    are now runtime parameters in the extern namelist as `eos_ttol`
    and `eos_dtol`.

## 18.08

  * the unit tests (`test_react`, `test_sdc`, and `test_eos`) have been
    ported from the Fortran to C++ build system in AMReX.  This will
    allow us to test the GPU framework in AMReX.

## 18.07

  * added CUDA support to the VODE90 integrator, the helmeos, and the
    networks `aprox13`, `aprox19`, `aprox21`, `ignition_simple`,
    `C-burn-simple`, `URCA-simple`.

  * Ported the unit test frameworks to FBoxLib

## 18.05

  * lots of documentation updates

  * some fixes to the numerical Jacobian involving X vs. Y (#100)

  * a new `subCh` network for He burning was added.

  * implemented the new c12(a,g)o16 nuclear reaction rate and its
    corresponding inverse from the work of Deboer et al. 2017 (Rev Mod
    Phys 89, 035007, 2017).  To use the new rate, user must set
    `use_c12ag_deboer17` to `true`.  This rate is only usable in the
    `aprox13`, `aprox19`, `aprox21`, and `iso7` reaction rate
    networks. (#56)

  * a routine `util/cj_detonation` was added to compute the
     Chapman-Jouguet detonation velocity for any of the networks

  * the burn retry strategy now sticks with the current integrator and
    uses looser tolerances before switching to a different integrator. (#96)

## 18.04

   * pynucastro (https://github.com/pynucastro/pynucastro) can now
     generate reaction networks compatible with StarKiller.  See the
     `subch` network.

## 17.11

   * a new option to boost the reaction rates has been added
     to the integrators (#64)

   * we now disable some composition derivatives in the EOS
     by default, for performance and memory reasons.  They can
     be re-enabled by defining the preprocessor variable
     `EXTRA_THERMO` (#59)

## 17.10

  * the compositional derivatives are no longer available
    by default from the EOS.  To get these, set the preprocessor
    variable `EXTRA_THERMO`.  This change was done for performance
    reasons.

  * the `aprox19` and `aprox21` networks no longer use a numerical
    Jacobian by default, as this was found to result in some
    bad numerical issues in VODE (#49)

  * the maximum temperature for reactions, `MAX_TEMP`, is now
    an adjustable input parameter rather than being hardcoded
    at `1.d11`.

  * the Helmholtz EOS table is now read by the IO processor and
    broadcast to other processors (#53)


  * the VODE integrator now does some additional checks on the
    state to ensure consistency (#47)

## 17.09

  * a new rety mechanism was implemented that allows a different
    integrator to be used if the primary integrator fails (#39)

  * the electron Ni56 electron capture rates and energy losses
    were updated from Mazurek (1973) to LMP (2000).  Thanks to
    Carl Fields for this contribution.  (#40)

## 17.08

  * fix to `aprox21` from Aron Michel (HITS) that fills in missing
    reactions

  * updated the helmholtz EOS to use the latest table from Frank
    Timmes (in particular, this is now denser, with 2x points in T and
    rho dimensions).  If you copied the old table, you need to make sure
    you are using the new table now.

  * add stellar conductivities from Frank Timmes

## 17.06

  * a new Fortran 90 port of VODE has been added

  * the unit tests now require AMReX instead of BoxLib to build

## 17.01

  * we've removed the option to integrate molar fractions and instead
    the ODE system always operates on mass fractions (the networks
    return derivatives of molar fractions and they are automatically
    converted).

## 16.12

  * a new unit test, `test_sdc`, was created to test the SDC interface
    to the networks

  * we now rely on the network module to provide `aion_inv` (1/aion)

  * the VODE integrator now supports SDC integration

## 16.09

  * num_rate_groups is now a property of the individual networks

  * a new integration method, Rosenbrock, was added to the BS
    option (set `ode_method`)

  * the number of RHS and Jac evaluations is now passed out
    of the burner through the burn_t type for diagnostic and
    load-balancing use

  * support for spectral deferred correction coupling of the
    burner and hydro was added to the BS integrator

## 16.08

  * `Microphysics/eos/` has been renamed `Microphysics/EOS/` to better
    conform to the conventions used in Castro and Maestro

  * the User's Guide has been extensively updated

  * OpenMP and OpenACC directives have been added to the unit tests

  * the BS integrator's type, `bs_t`, has now contains a `burn_t`
    internally, simplifying the conversion from `bs_t` to call the
    `actual_rhs/jac`

  * the `rates()` component of `burn_t` was removed.  We no longer
    rely on rate caching

  * we now store the simulation time at the start of the burn as
    t0 in the rpar storage to use as an offset to the integration
    time

  * the species derivatives (dh/dX and de/dX) and enthalpy were
    removed from the `burn_t`

  * a new option to integrate of X instead of Y was added
    (`integrate_molar_fraction = F`)

  * integration of networks with `nspec_evolve` < `nspec` were fixed
    to now apply the algrebic constraint relating mass fractions
    through a new `update_unevolved_species()` function

  * the electron capture rate on Ni56 used by `aprox19` and `aprox21` was
    fixed

  * the BS integrator can now use the initial timestep estimation
    algorithm that VODE uses (`use_timestep_estimator = T`)

  * a centered difference numerical Jacobian option was added


## 16.07

  * we now use `MICROPHYSICS_HOME` instead of `MICROPHYSICS_DIR` as the
    environment variable to point to the `Microphysics/` directory.

  * there are now two standalone unit tests, `test_react` and `test_eos`
    that don't need Maestro or Castro to compile.

  * new burn mode that limits numerically unstable burning.

  * `UsersGuide/` was renamed to `Docs/` to be consistent with the other
    BoxLib codes

  * the energy equation now uses an offset to help with the BS ODE
    integration convergence

  * the runtime parameter `small_x` now is owned by the network
