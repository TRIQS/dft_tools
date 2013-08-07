.. index:: full charge self consistency

Full charge self consistency
============================

.. warning::
  Before using this tool, you should be familiar with the band-structure package :program:`Wien2k`, since
  the calculation is controlled by the :program:`Wien2k` scripts! See also the :download:`dmftproj tutorial<TutorialDmftproj.pdf>`.

In order to do charge self-consistent calculations, we have to tell the band structure program about the
changes in the charge density due to correlation effects. In the following, we discuss how to use the 
:program:`TRIQS` tools in combination with the :program:`Wien2k` program, although an extension to other 
codes is also possible.

We can use the DMFT script as introduced in sections :ref:`LDADMFTmain` and :ref:`advanced`, with a few simple 
modifications. First, in order to be compatible with the :program:`Wien2k` standards, the DMFT script has to be 
named ``case.py``, where `case` is the name of the :program:`Wien2k` calculation, see the section 
:ref:`interfacetowien` for details. Then we set the variable 
`lda_filename` dynamically::

  import os
  lda_filename = os.getcwd().rpartition('/')[2]

This sets the `lda_filename` to the name of the current directory. The reminder of the scripts is completely the 
same as in one-shot calculations. Only at the very end we have to calculate the modified charge density,
and store it in a format such that :program:`Wien2k` can read it. Therefore, after the DMFT loop that we saw in the 
previous section, we symmetrise the self energy, and recalculate the impurity Green function::

  SK.symm_deg_gf(S.Sigma,orb=0)
  S.G <<= inverse(S.G0) - S.Sigma
  S.G.invert()

These steps are not necessary, but can help to reduce fluctuation of the total energy. 
Now we calculate the modified charge density::

  # find exact chemical potential
  SK.put_Sigma(Sigma_imp = [ S.Sigma ])
  chemical_potential = SK.find_mu( precision = 0.000001 )
  dN,d = SK.calc_density_correction(filename = lda_filename+'.qdmft')
  SK.save()

First we find the chemical potential with high precision, and after that the routine 
``SK.calc_density_correction(filename)`` calculates the density matrix including correlation effects. The result
is stored in the file `Filename`, which is later read by the :program:`Wien2k` program. The last statement saves 
the chemical potential into the hdf5 archive.
We need also the correlation energy, which we evaluate by the Migdal formula::

  correnerg = 0.5 * (S.G * S.Sigma).total_density()

From this value, we have to substract the double counting energy::

  correnerg -= SK.dc_energ[0]

and save this value into the file::

  if (mpi.is_master_node()):
    f=open(lda_filename+'.qdmft','a')
    f.write("%.16f\n"%correnerg)
    f.close()

The above steps are valid for a calculation with only one correlated atom in the unit cell, the most likely case
where you will apply this method. That is the reason why we give the index `0` in the list `SK.dc_energ`.
If you have more than one correlated atom in the unit cell, but all of them
are equivalent atoms, you have to multiply the `correnerg` by their multiplicity, before writing it to the file.
The multiplicity is easily found in the main input file of the :program:`Wien2k` package, i.e. `case.struct`. In case of
non-equivalent atoms, the correlation energy has to be calculated for all of them separately (FOR EXPERTS ONLY).

As mentioned above, the calculation is controlled by the :program:`Wien2k` scripts and not by :program:`python` 
routines. Therefore, you start your calculation for instance by::

  me@home $ run -qdmft -i 10

The flag `-qdmft` tells the script, that the density matrix including correlation effects is read from the `case.qdmft`
file, and 10 self-consitency iterations are done. If you run the code on a parallel machine, you can specify the number of 
nodes that are used::

  me@home $ run -qdmft -np 64 -i 10

with the `-np` flag. In that case, you have to give the proper `MPI` execution statement, e.g. `mpiexec`, in the `run_lapw` script, 
see the corresponding :program:`Wien2k` documentation. In many cases it is advisable to start from a converged one-shot 
calculation.

For practical purposes, you keep the number of DMFT loops within one DFT cycle low, or even to `loops=1`. If you encouter 
unstable convergence, you have to adjust the parameters such as
`loops`, `mix`, or `Delta_mix` to improve the convergence.

In the next section, :ref:`LDADMFTtutorial`, we will see in a detailed
example, how such a self consistent calculation is performed.

