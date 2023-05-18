The following equations of state are provided:

* `gamma_law`: a general constant-gamma ideal gas.  You can
  set whether or not the gas is ionized via a runtime parameter.

* `helmholtz`: a general stellar equation of state consisting of
  ions, degenerate/relativitic electrons, and radiation.  This is Frank
  Timmes's EOS with a modified interface and thread-safety.

* `multigamma`: a gamma-law EOS where each species can have its own
  gamma (ratio of specific heats).

* `polytrope`: a simple polytropic EOS: `P = K rho**gamma`

* `stellarcollapse`: the nuclear equation of state from
  stellarcollapse.org

* `ztwd`: a zero-temperature white dwarf equation of state

* `primordial_chem`: a version of multigamma for primordial ISM chemistry
