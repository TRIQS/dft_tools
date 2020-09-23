.. highlight:: python

.. _singleshot:

Single-shot DFT+DMFT
====================

After having set up the hdf5 archive, we can now proceed to our first DFT+DMFT calculation.
It consists of initialization steps, and the actual DMFT self-consistency loop,
With the code snippets below you can build your own script and target
it to your needs. Little examples on :ref:`mixing <mixing>` and on
:ref:`restarting from a previous calculation <restartcalc>` at the end of this page
should also demonstrate how simple you can modify your own DMFT script. A full working
calculation for SrVO3 is discussed in the :ref:`next section <SrVO3>`.


Initialization of the calculation
---------------------------------

Before doing the actual calculation, we have to initialize all needed objects.
The first thing is the :class:`SumkDFT <dft.sumk_dft.SumkDFT>` class.
It contains all basic routines that are necessary to perform a summation in k-space
to get the local quantities used in DMFT. It is initialized by::

    from triqs_dft_tools.sumk_dft import *
    SK = SumkDFT(hdf_file = filename + '.h5')


Setting up the impurity solver
------------------------------

The next step is to setup an impurity solver. There are different
solvers available within the :ref:`TRIQS <triqslibs:welcome>` framework.
E.g. for :ref:`SrVO3 <SrVO3>`, we will use the hybridization
expansion :ref:`CTHYB solver <triqscthyb:welcome>`. Later on, we will
see also the example of the `Hubbard-I solver <https://triqs.github.io/triqs/1.4/applications/hubbardI/>`_.
They all have in common, that they are called by an uniform command::

    S.solve(params)

where :emphasis:`params` are the solver parameters and depend on the actual
solver. Setting up the :ref:`CTHYB solver <triqscthyb:welcome>` for SrVO3 is
discussed on the :ref:`next page <SrVO3>`. Here, let us now perform the DMFT
loop using the methods of :program:`DFTTools`, assuming that we have already
set up a working solver instance.


Doing the DMFT loop
-------------------

Having initialized the :class:`Sumk class <dft.sumk_dft.SumkDFT>`
and the solver, we can proceed with the actual DMFT part of the calculation.
We set up the loop over DMFT iterations and the self-consistency condition::

    n_loops = 15
    for iteration_number in range(n_loops) :        # start the DMFT loop
        SK.set_Sigma([ S.Sigma ])                   # Put self energy to the SumK class
        chemical_potential = SK.calc_mu()           # calculate the chemical potential for the given density
        S.G_iw << SK.extract_G_loc()[0]                     # extract the local Green function
        S.G0_iw << inverse(S.Sigma_iw + inverse(S.G_iw))    # finally get G0, the input for the solver

        S.solve(h_int=h_int, **p)                   # now solve the impurity problem

        dm = S.G_iw.density()                                           # Density matrix of the impurity problem
        SK.calc_dc(dm, U_interact=U, J_hund=J, orb=0, use_dc_formula=1) # Set the double counting term
        SK.save(['chemical_potential','dc_imp','dc_energ'])             # Save data in the hdf5 archive

These steps are enough for a basic DMFT Loop.
After the self-consistency steps, which lead to a new :math:`G^0(i\omega)`,
the impurity solver is called. Different to model calculations, we have to do a few
more steps after this, because of the double-counting correction. We first
calculate the density of the impurity problem. Then, the routine :meth:`calc_dc <dft.sumk_dft.SumkDFT.calc_dc>`
takes as parameters this density matrix, the Coulomb interaction, Hund's rule
coupling, and the type of double-counting that should be used. Possible values
for :emphasis:`use_dc_formula` are:

    * `0`: Full-localised limit (FLL)
    * `1`: DC formula as given in K. Held, Adv. Phys. 56, 829 (2007).
    * `2`: Around-mean-field (AMF)

At the end of the calculation, we can save the Green function and self energy into a file::

    from h5 import HDFArchive
    import triqs.utility.mpi as mpi
    if mpi.is_master_node():
        ar = HDFArchive("YourDFTDMFTcalculation.h5",'w')
        ar["G"] = S.G_iw
        ar["Sigma"] = S.Sigma_iw

These are the essential steps necessary for a one-shot DFT+DMFT calculation.
For a detailed description of the :class:`SumkDFT <dft.sumk_dft.SumkDFT>`
routines, see the :ref:`reference manual <reference>`. To perform full charge self-consistent calculations, there
are some more things to consider, which we will see :ref:`later on <full_charge_selfcons>`.

.. _restartcalc:


Restarting a calculation
------------------------

Often only a few DMFT iterations are performed first, and thus, it is desirable to
carry out further iterations, e.g. to improve on the convergence. With a little modification
at the initialization stage (before the DMFT loop) it is possible to detect if previous runs
are present, or if the calculation should start from scratch::

    previous_runs = 0
    previous_present = False
    if mpi.is_master_node():
        with HDFArchive(dft_filename+'.h5','a') as f:
            if 'dmft_output' in f:
                ar = f['dmft_output']
                if 'iterations' in ar:
                    previous_present = True
                    previous_runs = ar['iterations']
            else:
                f.create_group('dmft_output')
    
    previous_runs = mpi.bcast(previous_runs)
    previous_present = mpi.bcast(previous_present)


You can see from this code snippet, that removing the subgroup :emphasis:`dmft_results` from the
hdf file has the effect of reseting the calculation to the starting point. If there are previous
runs stored in the hdf5 archive, we can now load the self energy, the chemical potential and
double counting values of the last iteration::

    if previous_present:
        if mpi.is_master_node():
            with HDFArchive(dft_filename+'.h5','r') as ar:
                S.Sigma_iw << ar['dmft_output']['Sigma_iw']

        S.Sigma_iw << mpi.bcast(S.Sigma_iw)
        chemical_potential,dc_imp,dc_energ = SK.load(['chemical_potential','dc_imp','dc_energ'])
        SK.set_mu(chemical_potential)
        SK.set_dc(dc_imp,dc_energ)

The data is loaded only on the master node, and therefore we broadcast it to the slave nodes.
Be careful when storing the :emphasis:`iteration_number` as we also have to add the previous
iteration count::

    ar['dmft_output']['iterations'] = iteration_number + previous_runs

.. _mixing:


Mixing
------

In some cases a mixing of two consecutive self energies (or alternatively two hybridization
functions) can be necessary in order to ensure convergence::

    mix = 0.8 # mixing factor
    if (iteration_number>1 or previous_present):
        if mpi.is_master_node():
            with HDFArchive(dft_filename+'.h5','r') as ar:
                mpi.report("Mixing Sigma and G with factor %s"%mix)
                S.Sigma_iw << mix * S.Sigma_iw + (1.0-mix) * ar['dmft_output']['Sigma_iw']
                S.G_iw << mix * S.G_iw + (1.0-mix) * ar['dmft_output']['G_iw']
        S.G_iw << mpi.bcast(S.G_iw)
        S.Sigma_iw << mpi.bcast(S.Sigma_iw)

In this little piece of code, which should be placed after calling the solver, two consecutive
self energies are linearly mixed with the factor :emphasis:`mix`. Of course, it is possible
to implement more advanced mixing schemes (e.g. Broyden's methods), however, in most cases
simple linear mixing or even no mixing is sufficient for a reasonably fast convergence.
