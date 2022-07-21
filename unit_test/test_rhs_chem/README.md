# `test_rhs`

This is a unit test that sets up a cube of data (rho, T, and X varying
along dimensions) and calls the RHS on it.

E.g. testing iso7 ...

```
~/d/M/u/test_rhs (iso7-cpp *$%>) ./main3d.gnu.ex inputs_iso7 do_cxx=0 amrex.fpe_trap_invalid=1
AMReX (20.04-59-g5ba96e09f85e) initialized
reading extern runtime parameters ...
 
  Initializing Helmholtz EOS and using Coulomb corrections.
   
   Run time = 0.046359805
   [The  Pinned Arena] space allocated (MB): 8
   [The  Pinned Arena] space used      (MB): 0
   AMReX (20.04-59-g5ba96e09f85e) finalized

~/d/M/u/test_rhs (iso7-cpp *$%>) ./main3d.gnu.ex inputs_iso7 do_cxx=1 amrex.fpe_trap_invalid=1
AMReX (20.04-59-g5ba96e09f85e) initialized
reading extern runtime parameters ...

 Initializing Helmholtz EOS and using Coulomb corrections.
  
  Run time = 0.03039193
  [The  Pinned Arena] space allocated (MB): 8
  [The  Pinned Arena] space used      (MB): 0
  AMReX (20.04-59-g5ba96e09f85e) finalized

~/d/M/u/test_rhs (iso7-cpp *$%>) fcompare.gnu.ex react_iso7_test_rhs.VODE.cxx/ react_iso7_test_rhs.VODE.fortran/

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
 Tdot                                          2005           9.795448465e-14
 Edot                               5.214358733e+10           7.830203405e-14
 J_helium-4_helium-4                2.793967724e-09           5.726429519e-16
 J_carbon-12_helium-4               9.313225746e-10           5.726450599e-16
 J_oxygen-16_helium-4               1.776356839e-15           1.849916399e-17
 J_neon-20_helium-4                 4.263256415e-14           4.608758669e-16
 J_magnesium-24_helium-4            4.263256415e-14           1.962564252e-15
 J_silicon-28_helium-4               3.60822483e-16           1.784402189e-15
 J_nickel-56_helium-4                             0                         0
 J_T_helium-4                                 52096           1.484699234e-13
 J_E_helium-4                       1.667521053e+12           1.460722174e-13
 J_helium-4_carbon-12               3.885780586e-16           1.381797231e-16
 J_carbon-12_carbon-12              3.885780586e-16           1.381797231e-16
 J_oxygen-16_carbon-12              3.885780586e-16           1.381797231e-16
 J_neon-20_carbon-12                1.447566072e-24           3.496975052e-15
 J_magnesium-24_carbon-12                         0                         0
 J_silicon-28_carbon-12                           0                         0
 J_nickel-56_carbon-12                            0                         0
 J_T_carbon-12                       0.007617950439           1.275406935e-14
 J_E_carbon-12                               198784           1.022949625e-14
 J_helium-4_oxygen-16               5.551115123e-16           1.443234513e-18
 J_carbon-12_oxygen-16                            0                         0
 J_oxygen-16_oxygen-16              5.551115123e-16           1.443234513e-18
 J_neon-20_oxygen-16                5.551115123e-16           1.443234513e-18
 J_magnesium-24_oxygen-16                         0                         0
 J_silicon-28_oxygen-16             7.623063646e-43           5.046983552e-16
 J_nickel-56_oxygen-16                            0                         0
 J_T_oxygen-16                         0.0793762207            1.47088244e-15
 J_E_oxygen-16                              1590784           9.060683304e-16
 J_helium-4_neon-20                 9.094947018e-13           1.777795711e-15
 J_carbon-12_neon-20                              0                         0
 J_oxygen-16_neon-20                              0                         0
 J_neon-20_neon-20                  9.094947018e-13           1.777795711e-15
 J_magnesium-24_neon-20             9.094947018e-13           1.777795711e-15
 J_silicon-28_neon-20                             0                         0
 J_nickel-56_neon-20                              0                         0
 J_T_neon-20                              35.859375           2.537611348e-13
 J_E_neon-20                              814219264           1.771030155e-13
 J_helium-4_magnesium-24            9.769962617e-15           1.725575743e-15
 J_carbon-12_magnesium-24                         0                         0
 J_oxygen-16_magnesium-24                         0                         0
 J_neon-20_magnesium-24                           0                         0
 J_magnesium-24_magnesium-24           9.769962617e-15           1.725575743e-15
 J_silicon-28_magnesium-24           9.769962617e-15           1.725575743e-15
 J_nickel-56_magnesium-24                         0                         0
 J_T_magnesium-24                       0.869140625           5.184470245e-13
 J_E_magnesium-24                          25452544           4.666676799e-13
 J_helium-4_silicon-28                            0                         0
 J_carbon-12_silicon-28                           0                         0
 J_oxygen-16_silicon-28                           0                         0
 J_neon-20_silicon-28                             0                         0
 J_magnesium-24_silicon-28                         0                         0
 J_silicon-28_silicon-28                          0                         0
 J_nickel-56_silicon-28                           0                         0
 J_T_silicon-28                     2.220446049e-15           5.460341819e-16
 J_E_silicon-28                                   0                         0
 J_helium-4_nickel-56                             0                         0
 J_carbon-12_nickel-56                            0                         0
 J_oxygen-16_nickel-56                            0                         0
 J_neon-20_nickel-56                              0                         0
 J_magnesium-24_nickel-56                         0                         0
 J_silicon-28_nickel-56                           0                         0
 J_nickel-56_nickel-56                            0                         0
 J_T_nickel-56                      4.440892099e-15           5.457339125e-16
 J_E_nickel-56                                    0                         0
 J_helium-4_T                       1.517883041e-18           8.597194144e-16
 J_carbon-12_T                      5.421010862e-19           9.211798353e-16
 J_oxygen-16_T                      8.271806126e-25           1.031835673e-17
 J_neon-20_T                         5.95570041e-23           7.658531959e-16
 J_magnesium-24_T                    5.95570041e-23           1.557721422e-15
 J_silicon-28_T                     9.305781891e-25           1.682232215e-15
 J_nickel-56_T                                    0                         0
 J_T_T                              2.425909042e-05            1.91044887e-13
 J_E_T                                        777.5           1.882016793e-13
 J_helium-4_E                                     0                         0
 J_carbon-12_E                                    0                         0
 J_oxygen-16_E                                    0                         0
 J_neon-20_E                                      0                         0
 J_magnesium-24_E                                 0                         0
 J_silicon-28_E                                   0                         0
 J_nickel-56_E                                    0                         0
 J_T_E                                            0                         0
 J_E_E                                            0                         0
 ```
