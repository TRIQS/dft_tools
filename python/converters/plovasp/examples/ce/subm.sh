#!/bin/bash -l
#
#SBATCH --job-name="JOBNAME"
#SBATCH --time=36:00:00
##SBATCH --nodes=8
#SBATCH --ntasks=8
#SBATCH --partition=qmac
##SBATCH --mem-per-cpu=1024
#SBATCH --output=job.%j.o
#SBATCH --error=job.%j.e

#======START=====
# source $HOME/.profile
# export KMP_AFFINITY=none

# module load slurm
source /cm/shared/apps/intel/composer_xe/composer_xe_2015.0.090/mkl/bin/mklvars.sh intel64
# module switch PrgEnv-gnu PrgEnv-intel
# module load intel/14.0.2.144
#module load intel/mkl/64/11.2/2015.0.090

export OMP_NUM_THREADS=1

VASP_DIR=$HOME/Codes/vasp/build_20151202/vasp.5.4.2/build/std
#mpirun -np $SLURM_NTASKS $VASP_DIR/vasp
./sc_dmft.sh


