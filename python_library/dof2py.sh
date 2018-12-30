#!/usr/bin/sh
cd tmp_build_dir/o/3d.gnu.EXE
f2py -c -m _StarKillerMicrophysics *.o f90wrap_*.f90 -lstdc++
cd ../../..
mv tmp_build_dir/o/3d.gnu.EXE/*.so .
