# CVODE test react example

# Building on Groot

## Ignition Simple Network and CUDA
make -j 4 CUDA_ARCH=60 COMPILE_CUDA_PATH=/usr/local/cuda-9.2 USE_CUDA=TRUE COMP=PGI EOS_DIR=helmholtz INTEGRATOR_DIR=CVODE NETWORK_DIR=ignition_simple USE_MPI=FALSE USE_OMP=FALSE USE_ACC=FALSE CVODE_HOME=/home/dwillcox/dev-ode/cvode/instdir

## Ignition Simple Network and Serial
make -j 4 USE_CUDA=FALSE COMP=GNU EOS_DIR=helmholtz INTEGRATOR_DIR=CVODE NETWORK_DIR=ignition_simple USE_MPI=FALSE USE_OMP=FALSE USE_ACC=FALSE CVODE_HOME=/home/dwillcox/dev-ode/cvode/instdir USE_CUDA_CVODE=FALSE

