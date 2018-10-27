#!/usr/bin/sh
f90wrap -k .f90wrap_kind_map -m StarKillerMicrophysics tmp_build_dir/f/3d.gnu.EXE/F90PP-actual_burner.F90 tmp_build_dir/f/3d.gnu.EXE/F90PP-actual_network.F90 tmp_build_dir/f/3d.gnu.EXE/F90PP-burn_type.F90 tmp_build_dir/f/3d.gnu.EXE/F90PP-eos_type.F90 tmp_build_dir/f/3d.gnu.EXE/F90PP-eos.F90 tmp_build_dir/f/3d.gnu.EXE/F90PP-network.F90 tmp_build_dir/f/3d.gnu.EXE/F90PP-microphysics.F90 tmp_build_dir/f/3d.gnu.EXE/F90PP-extern.F90 starkiller_initialization.f90
mv f90wrap* tmp_build_dir/o/3d.gnu.EXE/.

