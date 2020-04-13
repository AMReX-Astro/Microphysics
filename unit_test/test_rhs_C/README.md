# `test_rhs`

This is a unit test that sets up a cube of data (rho, T, and X varying
along dimensions) and calls the RHS on it.

E.g. testing iso7 ...

```
~/d/M/u/test_rhs_C (iso7-cpp *+$%=) ./main3d.gnu.ex inputs_iso7 do_cxx=0 amrex.fpe_trap_invalid=1
AMReX (20.04-59-g5ba96e09f85e) initialized
reading extern runtime parameters ...
 
 Initializing Helmholtz EOS and using Coulomb corrections.
 
Run time = 0.014140926
[The  Pinned Arena] space allocated (MB): 8
[The  Pinned Arena] space used      (MB): 0
AMReX (20.04-59-g5ba96e09f85e) finalized

~/d/M/u/test_rhs_C (iso7-cpp *+$%=) ./main3d.gnu.ex inputs_iso7 do_cxx=1 amrex.fpe_trap_invalid=1
AMReX (20.04-59-g5ba96e09f85e) initialized
reading extern runtime parameters ...
 
 Initializing Helmholtz EOS and using Coulomb corrections.
 
Run time = 0.011533634
[The  Pinned Arena] space allocated (MB): 8
[The  Pinned Arena] space used      (MB): 0
AMReX (20.04-59-g5ba96e09f85e) finalized

~/d/M/u/test_rhs_C (iso7-cpp *+$%=) fcompare.gnu.ex react_iso7_test_rhs.VODE.cxx react_iso7_test_rhs.VODE.fortran/

            variable name            absolute error            relative error
                                        (||A - B||)         (||A - B||/||A||)
 ----------------------------------------------------------------------------
 level = 0
 density                                          0                         0
 temperature                                      0                         0
 Ydot_helium-4                      1.746229827e-10           6.135416039e-16
 Ydot_carbon-12                     5.820766091e-11           6.135483795e-16
 Ydot_oxygen-16                     3.469446952e-18           1.445247187e-18
 Ydot_neon-20                       2.664535259e-15           1.152189667e-15
 Ydot_magnesium-24                  2.664535259e-15           1.869108811e-15
 Ydot_silicon-28                    2.428612866e-17           1.830156091e-15
 Ydot_nickel-56                                   0                         0
 Xold_helium-4                                    0                         0
 Xold_carbon-12                                   0                         0
 Xold_oxygen-16                                   0                         0
 Xold_neon-20                                     0                         0
 Xold_magnesium-24                                0                         0
 Xold_silicon-28                                  0                         0
 Xold_nickel-56                                   0                         0
 Tdot                                    1714760544            8.37748127e-08
 Edot                               5.578806411e+16            8.37748127e-08
```
