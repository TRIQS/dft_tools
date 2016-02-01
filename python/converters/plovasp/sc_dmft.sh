#/bin/bash

NPROC=4

rm -f vasp.lock
stdbuf -o 0 mpirun -np $NPROC ~/Codes/vasp/build_20151202/vasp.5.4.2/build/std/vasp &

mpirun -np $NPROC ./run_build.sh sc_dmft.py $(jobs -p) || kill %1

