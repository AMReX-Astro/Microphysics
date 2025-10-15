The following equations of state are provided:

* `breakout`: this is essentially the same as `gamma_law`, but it sets
  the mean molecular weight based on ``eos_t aux[]`` instead of the
  mass fractions.

* `gamma_law` : a general constant-gamma ideal gas.  You can set
  whether or not the gas is ionized via a runtime parameter.

* `helmholtz`: a general stellar equation of state consisting of ions,
  degenerate/relativitic electrons, and radiation.  This is Frank
  Timmes's EOS with a modified interface and thread-safety.

* `metal_chem`: a multi-gamma EOS for metal ISM chemistry.

* `multigamma`: a gamma-law EOS where each species can have its own
  gamma (ratio of specific heats).

* `polytrope`: a simple polytropic EOS: `P = K rho**gamma`

* `primordial_chem`: a version of multigamma for primordial ISM
  chemistry

* `rad_power_law`: an artificial EOS for radiation transport test
  problems.

* `tillotson`: an EOS for hypevelocity impacts based on Tillotson
  (1962)

* `ztwd`: a zero-temperature white dwarf equation of state

