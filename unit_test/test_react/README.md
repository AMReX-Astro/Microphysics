# `test_react`

This is a unit test that sets up a cube of data (rho, T, and X varying
along dimensions) and calls the burner on it.  You can specify the integrator
via `INTEGRATOR_DIR` and the network via `NETWORK_DIR` in the `GNUmakefile`

## Status

| network | VODE | VBDF | BS | Rosenbrock | VBDF-ACC | BS-ACC |
|---------|------|------|----|------------|----------|--------|
| aprox13 | works<br>(6:673:5173) | | | | | |

