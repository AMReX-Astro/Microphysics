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


## Running on GPUs with VODE90

To run a GPU test with the VODE90 integrator and aprox13, do:

```
make -j COMP=PGI USE_CUDA=TRUE AMREX_USE_CUDA=TRUE USE_GPU_PRAGMA=TRUE NETWORK_DIR=aprox13 EOS_DIR=helmholtz
```

Tested with PGI 18.10 and CUDA 9.2.148.

To run in a different directory, the following files are necessary to
run in addition to the executable:

- `helm_table.dat` (will be symlinked in at compile time)
- `inputs_aprox13`
- `probin.aprox13`
- `xin.aprox13`

Then in an interactive session, (e.g. on Summit), do:

```
jsrun -n 1 -a 1 -g 1 ./[executable] inputs_aprox13
```
