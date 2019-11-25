#!/usr/bin/sh
# Usage in an AMReX makefile (after including Make.defs):
#       $ sh dof2py.sh $(objEXETempDir)
# First argument: AMReX temporary build directory holding compiled object files
cd $1
f2py3 -c -m _StarKillerMicrophysics *.o f90wrap_*.f90 -lstdc++
cd ../../..
mv $1/*.so .
