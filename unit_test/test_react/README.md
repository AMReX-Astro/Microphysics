# `test_react`

This is a unit test that sets up a cube of data (rho, T, and X varying
along dimensions) and calls the burner on it.  You can specify the integrator
via `INTEGRATOR_DIR` and the network via `NETWORK_DIR` in the `GNUmakefile`

## CPU Status

This table summarizes tests run with gfortran.


| network | VODE                  | VBDF                                 | BS                      | Rosenbrock |
|---------|-----------------------|--------------------------------------|-------------------------|------------|
| aprox13 | works<br>(6:673:5173) | crashes unless<br>min density raised | works<br>(126:858:4678) | crashes    |





## OpenACC Status
