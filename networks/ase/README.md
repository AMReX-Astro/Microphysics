# `ase`

This network is specifically tailored for self-consistent NSE evolution.
This network achieves machine precision accuracy when comparing equilibrium
composition from direct integration and NSE composition. 
This can be demonstrated using the `nse_compatibility` script.

This is similar to `he-burn-19a` except:

* N14 is not included.

* Every forward rate has a corresponding reverse rate calculated from
  detailed balance. This includes the reverse rates for C12+C12,
  C12+O16, and O16+O16 as they're removed in other `he-burn` nets.

* Q-value is recomputed in `DerivedRate` for better consistency
  with NSE.
