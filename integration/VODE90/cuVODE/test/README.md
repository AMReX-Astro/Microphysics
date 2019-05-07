# Using cuVODE

This example program shows how to use cuVODE in an AMReX application.

To build the example, you should have the environment variable
`AMREX_HOME` defined and then do `make`.

To use cuVODE with your own problem, you will need to write your own
versions of `cuvode_parameters.F90`, `rpar.F90`, `vode_rhs.F90`, and
`react_zones.F90`. The first three are expected by cuVODE and define
parameters of your ODE system. The last is what calls the `dvode`
integrator and will probably look similar to the ODE driver routine in
your application code that integrates ODEs in a loop over zones.

# Building with CUDA

To build the example with CUDA Fortran, use the PGI compiler and do:

```
make -j COMP=PGI USE_CUDA=TRUE AMREX_USE_CUDA=TRUE USE_GPU_PRAGMA=TRUE CUDA_VERSION=9.0
```

This was tested with PGI 18.10 and CUDA 9.2.148.

