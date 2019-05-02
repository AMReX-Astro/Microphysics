# CVODE test react example

This test packs the ODEs corresponding to the reaction networks in
individual cells in FABs together into a single vector of ODEs to
integrate with CVODE.

To use batched QR provided by CUDA cuSOLVER, build CVODE from the
master branch here:
https://github.com/dwillcox/cvode-4.0.0-development

In the build lines below, `CVODE_HOME` is the CVODE installation
directory, containing the `include` and `lib` subdirectories.

# Building on Summit

## Aprox13 + CVODE serial

The following compiles in the interface to serial CVODE in `react_serial.cpp`:

```
make -j USE_CUDA=FALSE COMP=GNU EOS_DIR=helmholtz INTEGRATOR_DIR=CVODE NETWORK_DIR=aprox13 USE_MPI=FALSE USE_OMP=FALSE USE_ACC=FALSE CVODE_HOME=/ccs/home/dwillcox/run-cuda-vode-cpp/cvode-cusolver/instdir USE_CVODE_CUSOLVER=FALSE USE_GPU_PRAGMA=FALSE USE_SPARSE_STOP_ON_OOB=FALSE DEBUG=FALSE USE_CUDA_CVODE=FALSE
```

## Aprox13 + CUDA + CVODE SPGMR solver

The following compiles in the interface to CUDA CVODE with the SPGMR linear solver in `react_cuda.cpp`:

```
make -j USE_CUDA=TRUE COMP=PGI EOS_DIR=helmholtz INTEGRATOR_DIR=CVODE NETWORK_DIR=aprox13 USE_MPI=FALSE USE_OMP=FALSE USE_ACC=FALSE CVODE_HOME=/ccs/home/dwillcox/run-cuda-vode-cpp/cvode-cusolver/instdir USE_CVODE_CUSOLVER=FALSE USE_GPU_PRAGMA=TRUE USE_SPARSE_STOP_ON_OOB=FALSE DEBUG=FALSE
```

## Aprox13 + CUDA + CVODE with cuSOLVER

The following compiles in the interface to CUDA CVODE with the cuSOLVER batched QR linear solver in `react_cuda_cusolver.cpp`:

```
make -j COMP=PGI USE_MPI=FALSE USE_OMP=FALSE USE_CUDA=TRUE USE_CUDA_CVODE=TRUE USE_CVODE_CUSOLVER=TRUE AMREX_USE_CUDA=TRUE USE_GPU_PRAGMA=TRUE USE_SPARSE_STOP_ON_OOB=FALSE NETWORK_DIR=aprox13 CVODE_HOME=/ccs/home/dwillcox/run-cuda-vode-cpp/cvode-cusolver/instdir
```

# Building on Groot

## Ignition Simple Network and CUDA

```
make -j 4 CUDA_ARCH=60 COMPILE_CUDA_PATH=/usr/local/cuda-9.2 USE_CUDA=TRUE COMP=PGI EOS_DIR=helmholtz INTEGRATOR_DIR=CVODE NETWORK_DIR=ignition_simple USE_MPI=FALSE USE_OMP=FALSE USE_ACC=FALSE CVODE_HOME=/home/dwillcox/dev-ode/cvode/instdir
```

## Ignition Simple Network and Serial

```
make -j 4 USE_CUDA=FALSE COMP=GNU EOS_DIR=helmholtz INTEGRATOR_DIR=CVODE NETWORK_DIR=ignition_simple USE_MPI=FALSE USE_OMP=FALSE USE_ACC=FALSE CVODE_HOME=/home/dwillcox/dev-ode/cvode/instdir USE_CUDA_CVODE=FALSE
```

# Comparing with test_react

For GPUs, this can be compared with the VODE90 integrator in the
test_react test located in ../test_react

To build the test_react test, in ../test_react, compile as:

```
make -j COMP=PGI DEBUG=FALSE USE_MPI=FALSE USE_OMP=FALSE USE_ACC=FALSE USE_CUDA=TRUE NETWORK_DIR=aprox13 INTEGRATOR_DIR=VODE90 EOS_DIR=helmholtz CUDA_VERSION=9.0
```

To run in a different directory, for the test_react and
test_cvode_react cases, the following files are necessary to run in
addition to the executables:

- `helm_table.dat` (will be symlinked in at compile time)
- `inputs_aprox13`
- `probin.aprox13`
- `xin.aprox13`

To run either the test_react or test_cvode_react cases, do:

```
./[executable] inputs_aprox13
```

On Summit, use:

```
jsrun -n 1 -a 1 -g 1 ./[executable] inputs_aprox13
```


