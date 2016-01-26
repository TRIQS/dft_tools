.. highlight:: python

#######################
Charge self-consistency
#######################

Introduction
************

Running a DFT+DMFT calculation in the charge self-consistent mode
requires additional process management whereby the control is given
to the VASP and to the TRIQS code in an alternating manner.
It is implemented by a lock system. A VASP process is launched in
the background and a self-consistency (SC) script in the foreground.
Once VASP reaches the point where the projectors are generated 
it creates a lock file `vasp.lock` and waits until the lock file is
removed. The SC script, in turn, waits for the VASP process and once
the lock file is created it starts a DMFT iteration. The DMFT iteration
must finish by generating a Kohn-Sham (KS) density matrix (file `GAMMA`)
and removing the lock file. The VASP process then reads in `GAMMA`
and proceeds with the next iteration.

The DMFT iteration is performed using a user-defined script.
However, the control part is maintained by two universal scripts:
`sc_dmft.sh` and `sc_dmft.py`. 

The first script, `sc_dmft.sh`, launches the VASP process
in the background mode, takes its pid, and launches `sc_dmft.py`
in the foreground mode.
Both processes are run within an MPI environment with an appropriate
number of MPI nodes.

The second script, `sc_dmft.py`, is responsible for the charge self-consitency
logic described in the first paragraph. It also combines the total
energy contributions coming from the DFT and DMFT parts.
The script imports a user-defined DMFT script `dmft_cycle.py` which
must produce KS density-matrix file `GAMMA` for the next DFT iteration
and also must return the correlation (including the double counting) energy
of the impurity model as well as a correlation correction to the 
DFT band energy (resulting from the difference between the bare
and DMFT density matrices).


