#!/usr/bin/sh
# Usage in an AMReX makefile (after including Make.defs):
#       $ sh dowrap.sh $(f77EXETempDir) $(objEXETempDir)
# First argument: AMReX temporary build directory holding preprocessed Fortran source files
# Second argument: AMReX temporary build directory holding compiled object files
f90wrap -k .f90wrap_kind_map -m StarKillerMicrophysics $1/F90PP-actual_burner.F90 $1/F90PP-actual_network.F90 $1/F90PP-actual_rhs.F90 $1/F90PP-numerical_jacobian.F90 $1/F90PP-burn_type.F90 $1/F90PP-eos_type.F90 $1/F90PP-eos.F90 $1/F90PP-network.F90 $1/F90PP-microphysics.F90 $1/F90PP-extern.F90 $1/F90PP-sneut5.F90 $1/F90PP-screen.F90 $1/F90PP-stellar_conductivity.F90 starkiller_initialization.f90 $1/F90PP-integrator.F90
mv f90wrap* $2/.
