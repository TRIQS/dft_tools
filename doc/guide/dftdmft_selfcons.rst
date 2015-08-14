.. _full_charge_selfcons:

Full charge self consistency
============================

Wien2k + dmftproj
-----------------


.. warning::
  Before using this tool, you should be familiar with the band-structure package :program:`Wien2k`, since
  the calculation is controlled by the :program:`Wien2k` scripts! Be
  sure that you also understand how :program:`dmftproj` is used to
  construct the Wannier functions. For this step, see either sections
  :ref:`conversion`, or the extensive :download:`dmftproj manual<images_scripts/TutorialDmftproj.pdf>`.

In order to do charge self-consistent calculations, we have to tell the band structure program about the
changes in the charge density due to correlation effects. In the following, we discuss how to use the 
:program:`TRIQS` tools in combination with the :program:`Wien2k` program.

We can use the DMFT script as introduced in section :ref:`singleshot`,
with just a few simple 
modifications. First, in order to be compatible with the :program:`Wien2k` standards, the DMFT script has to be 
named :file:`case.py`, where `case` is the place holder name of the :program:`Wien2k` calculation, see the section 
:ref:`conversion` for details. We can then set the variable `dft_filename` dynamically::

  import os
  dft_filename = os.getcwd().rpartition('/')[2]

This sets the `dft_filename` to the name of the current directory. The
remaining part of the script is identical to 
that for one-shot calculations. Only at the very end we have to calculate the modified charge density,
and store it in a format such that :program:`Wien2k` can read it. Therefore, after the DMFT loop that we saw in the 
previous section, we symmetrise the self energy, and recalculate the impurity Green function::

  SK.symm_deg_gf(S.Sigma,orb=0)
  S.G_iw << inverse(S.G0_iw) - S.Sigma_iw
  S.G_iw.invert()

These steps are not necessary, but can help to reduce fluctuations in the total energy. 
Now we calculate the modified charge density::

  # find exact chemical potential
  SK.put_Sigma(Sigma_imp = [ S.Sigma_iw ])
  chemical_potential = SK.calc_mu( precision = 0.000001 )
  dN, d = SK.calc_density_correction(filename = dft_filename+'.qdmft')
  SK.save(['chemical_potential','dc_imp','dc_energ'])

First we find the chemical potential with high precision, and after that the routine 
``SK.calc_density_correction(filename)`` calculates the density matrix including correlation effects. The result
is stored in the file `dft_filename.qdmft`, which is later read by the :program:`Wien2k` program. The last statement saves 
the chemical potential into the hdf5 archive.

We need also the correlation energy, which we evaluate by the Migdal formula::

  correnerg = 0.5 * (S.G_iw * S.Sigma_iw).total_density()

Other ways of calculating the correlation energy are possible, for
instance a direct measurment of the expectation value of the
interacting hamiltonian. However, the Migdal formula works always,
independent of the solver that is used to solve the impurity problem.
From this value, we substract the double counting energy::

  correnerg -= SK.dc_energ[0]

and save this value in the file, too::

  if (mpi.is_master_node()):
    f=open(dft_filename+'.qdmft','a')
    f.write("%.16f\n"%correnerg)
    f.close()

The above steps are valid for a calculation with only one correlated atom in the unit cell, the most likely case
where you will apply this method. That is the reason why we give the index `0` in the list `SK.dc_energ`.
If you have more than one correlated atom in the unit cell, but all of them
are equivalent atoms, you have to multiply the `correnerg` by their multiplicity before writing it to the file.
The multiplicity is easily found in the main input file of the :program:`Wien2k` package, i.e. `case.struct`. In case of
non-equivalent atoms, the correlation energy has to be calculated for
all of them separately and summed up.

As mentioned above, the calculation is controlled by the :program:`Wien2k` scripts and not by :program:`python` 
routines. You should think of replacing the lapw2 part of the
:program:`Wien2k` self-consistency cycle by

  |   `lapw2 -almd`
  |   `dmftproj`
  |   `pytriqs case.py`
  |   `lapw2 -qdmft`

In other words, for the calculation of the density matrix in lapw2, we
add the DMFT corrections through our python scripts.
Therefore, at the command line, you start your calculation for instance by:

  `me@home $ run -qdmft 1 -i 10`

The flag `-qdmft` tells the :program:`Wien2k` script that the density
matrix including correlation effects is to be read in from the
`case.qdmft` file, and that you want the code to run on one computing
core only. Moreover, we ask for 10 self-consistency iterations are to be
done. If you run the code on a parallel machine, you can specify the
number of nodes to be used:

  `me@home $ run -qdmft 64 -i 10`

In that case, you will run on 64 computing cores. As standard setting,
we use `mpirun` as the proper MPI execution statement. If you happen
to have a differnet, non-standard MPI setup, you have to give the
proper MPI execution statement, in the `run_lapw` script (see the  
corresponding :program:`Wien2k` documentation).

In many cases it is advisable to start from a converged one-shot 
calculation. For practical purposes, you keep the number of DMFT loops
within one DFT cycle low, or even to `loops=1`. If you encouter
unstable convergence, you have to adjust the parameters such as
the number of DMFT loops, or some mixing of the self energy to improve
the convergence. 

In the section :ref:`DFTDMFTtutorial` we will see in a detailed
example how such a self-consistent calculation is performed from scratch.


Other DFT codes
---------------

The extension to other DFT codes is straight forward. As described
here, one needs to implement the correlated density matrix to be used
for the calculation of the charge density. This implementation will of
course depend on the DFT package, and might be easy to do or a quite
involved project. The formalism, however, is straight forward.
