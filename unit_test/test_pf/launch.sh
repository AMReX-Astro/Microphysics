#!/bin/bash
#SBATCH -N 1
#SBATCH -C haswell
#SBATCH -q debug
#SBATCH -t 00:00:15

#OpenMP settings:
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread


#run the application:
srun -n 1 -c 64 --cpu_bind=cores ./main3d.gnu.haswell.ex inputs_C_MG_burn
