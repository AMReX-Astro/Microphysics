# `test_react`

This is a unit test that sets up a cube of data (rho, T, and X varying
along dimensions) and calls the burner on it.  You can specify the integrator
via `INTEGRATOR_DIR` and the network via `NETWORK_DIR` in the `GNUmakefile`

## CPU Status

This table summarizes tests run with gfortran.


| network                | VODE                    | VODE90                  | VBDF                                 | BS                      | Rosenbrock    |
|------------------------|-------------------------|-------------------------|--------------------------------------|-------------------------|---------------|
| aprox13                | works<br>(6:1517:2393)  | works<br>(6:1517:23793) | crashes unless<br>min density raised | works<br>(126:858:4678) | crashes       |
| aprox19                | works<br>(27:152:8127)  | works<br>(27:152:8127)  | crashes                              | works<br>(72:297:2869)  | runs too long |
| triple_alpha_plus_cago | works<br>(6:32:1950)    | works<br>(6:32:1950)    | works<br>(4:37:1964)                 | works<br>(126:148:3044) | crashes       |
| ignition_chamulak      | works<br>(6:6:28)       | works<br>(6:6:28)       | works<br>(4:5:30)                    | works<br>(144:144:153)  | works (252:252:252) |




## OpenACC Status
