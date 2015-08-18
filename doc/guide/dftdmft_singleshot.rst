.. highlight:: python

.. _singleshot:

Single-shot DFT+DMFT
====================


After having set up the hdf5 archive, we can now do our DFT+DMFT calculation. It consists of
initialization steps, and the actual DMFT self-consistency loop, as is
discussed below. 

Initialisation of the calculation
---------------------------------

Before doing the calculation, we have to intialize all the objects that we will need. The first thing is the 
:class:`SumkDFT` class. It contains all basic routines that are necessary to perform a summation in k-space 
to get the local quantities used in DMFT. It is initialized by::

  from pytriqs.applications.dft.sumk_dft import *
  SK = SumkDFT(hdf_file = filename + '.h5')


Setting up the impurity solver
------------------------------

The next step is to setup an impurity solver. There are different
solvers available within the TRIQS framework. Below, we will discuss
the example of the hybridisation
expansion :ref:`CTHYB solver <triqscthyb:welcome>`. Later on, we will
see also the example of the Hubbard-I solver. They all have in common,
that they are called by a uniform command::
  
  S.solve(params)

where `params` are the solver parameters and depend on the actual
solver that is used. Before going into the details of the solver, let
us discuss in the next section how to perform the DMFT loop using
the methods of :program:`dft_tools`, assuming that we have set up a
working solver instance. 


Doing the DMFT loop
-------------------

Having initialized the SumK class and the Solver, we can proceed with the DMFT
loop itself. We have to set up the loop over DMFT
iterations and the self-consistency condition::

  n_loops = 5
  for iteration_number in range(n_loops) :            # start the DMFT loop

          SK.put_Sigma(Sigma_imp = [ S.Sigma ])              # Put self energy to the SumK class
          chemical_potential = SK.calc_mu()                  # calculate the chemical potential for the given density
          S.G_iw << SK.extract_G_loc()[0]                    # extract the local Green function
          S.G0_iw << inverse(S.Sigma_iw + inverse(S.G_iw))   # finally get G0, the input for the Solver

          S.solve(h_int=h_int, **p)                          # now solve the impurity problem

	  dm = S.G_iw.density()                                                 # Density matrix of the impurity problem  
          SK.calc_dc(dm, U_interact=U, J_hund=J, orb=0,	use_dc_formula=1)       # Set the double counting term
          SK.save(['chemical_potential','dc_imp','dc_energ'])                   # Save data in the hdf5 archive

These basic steps are enough to set up the basic DMFT Loop. For a detailed
description of the :class:`SumkDFT` routines, see the reference
manual.

After
the self-consistency steps (extracting a new :math:`G^0(i\omega)`),
the Anderson impurity problem is solved. 

Different to model calculations, we have to do a few
more steps after this, because of the double-counting correction. We first
calculate the density of the impurity problem. Then, the routine `calc_dc`
takes as parameters this density matrix, the Coulomb interaction, Hund's rule
coupling, and the type of double-counting that should be used. Possible values
for `use_dc_formula` are:

  * `0`: Full-localised limit
  * `1`: DC formula as given in K. Held, Adv. Phys. 56, 829 (2007).
  * `2`: Around-mean-field

At the end of the calculation, we can save the Greens function and self energy into a file::

  from pytriqs.archive import HDFArchive
  import pytriqs.utility.mpi as mpi
  if mpi.is_master_node():
      ar = HDFArchive("YourDFTDMFTcalculation.h5",'w')
      ar["G"] = S.G_iw
      ar["Sigma"] = S.Sigma_iw

This is it! 

These are the essential steps to do a one-shot DFT+DMFT calculation. 
For full charge-self consistent calculations, there are some more things 
to consider, which we will see later on.


A full DFT+DMFT calculation
---------------------------

We will discuss now how to set up a full working calculation,
including setting up the CTHYB solver, and specifying some more parameters
in order to make the calculation more efficient. Here, we
will see a more advanced example, which is also suited for parallel
execution. For the convenience of the user, we provide also two
working python scripts in this documentation. One for a calculation
using Kanamori definitions (:download:`dft_dmft_cthyb.py
<images_scripts/dft_dmft_cthyb.py>`) and one with a
rotational-invariant Slater interaction hamiltonian (:download:`dft_dmft_cthyb_slater.py
<images_scripts/dft_dmft_cthyb.py>`). The user has to adapt these
scripts to his own needs.

First, we load the necessary modules::

  from pytriqs.applications.dft.sumk_dft import *
  from pytriqs.gf.local import *
  from pytriqs.archive import HDFArchive
  from pytriqs.operators.util import *
  from pytriqs.applications.impurity_solvers.cthyb import *

The last two lines load the modules for the construction of the CTHYB
solver.

Then we define some parameters::

  dft_filename='SrVO3'
  U = 4.0
  J = 0.65
  beta = 40
  loops =  10                      # Number of DMFT sc-loops
  sigma_mix = 0.8                  # Mixing factor of Sigma after solution of the AIM
  dc_type = 1                      # DC type: 0 FLL, 1 Held, 2 AMF
  use_blocks = True                # use bloc structure from DFT input
  prec_mu = 0.0001

  # Solver parameters
  p = {}
  p["length_cycle"] = 200
  p["n_warmup_cycles"] = 2000
  p["n_cycles"] = 20000

Most of these parameters are self-explanatory. The first,
`dft_filename`, gives the filename of the input files. For more
details on the solver parameters, we refer the user to
the :ref:`CTHYB solver <triqscthyb:welcome>` documentation.

We assume that the conversion to the hdf5 archive is alreadz done. We
can check now in this archive, if previous runs are present, or if we have to start
from scratch::

  previous_runs = 0
  previous_present = False
  if mpi.is_master_node():
      f = HDFArchive(dft_filename+'.h5','a')
      if 'dmft_output' in f:
          ar = f['dmft_output']
          if 'iterations' in ar:
              previous_present = True
              previous_runs = ar['iterations']
      else:
          f.create_group('dmft_output')
      del f
  previous_runs    = mpi.bcast(previous_runs)
  previous_present = mpi.bcast(previous_present)


You can see in this code snipet, that all results of this calculation
will be stored in a separate subgroup in the hdf5 file, called
`dmft_output`. Removing this subgroup allows you to reset your
calculation to the starting point easily.

Now we can use all this information to initialise the :class:`SumkDFT` class::

  SK = SumkDFT(hdf_file=dft_filename+'.h5',use_dft_blocks=use_blocks)

The next step is to initialise the  :class:`Solver` class. It consist
of two steps

#. Calculating the multi-band interaction matrix, and setting up the
   interaction hamiltonian
#. Setting up the solver class

The first step is done using methods of
the :ref:`TRIQS <triqslibs:welcome>` library::

  n_orb = SK.corr_shells[0]['dim']
  l = SK.corr_shells[0]['l']
  spin_names = ["up","down"]
  orb_names = [i for i in range(n_orb)]
  # Use GF structure determined by DFT blocks:
  gf_struct = SK.gf_struct_solver[0]
  # Construct U matrix for density-density calculations:
  Umat, Upmat = U_matrix_kanamori(n_orb=n_orb, U_int=U, J_hund=J)

We assumed here that we want to use an interaction matrix with
Kanamori definitions of :math:`U` and :math:`J`. For
other choices (Slater interaction matrix for instance), and other
parameters, we refer to the reference manual 
of the :ref:`TRIQS <triqslibs:welcome>` library.

Next, we construct the hamiltonian and the solver::
  
  h_int = h_int_density(spin_names, orb_names, map_operator_structure=SK.sumk_to_solver[0], U=Umat, Uprime=Upmat)
  S = Solver(beta=beta, gf_struct=gf_struct)

As you see, we take only density-density interactions into
account. Other choices for the hamiltonian are

* h_int_kanamori
* h_int_slater

These two include full rotational invariant interactions. Again,
options can be found in the :ref:`TRIQS <triqslibs:welcome>` library
reference manual.


If there are previous runs stored in the hdf5 archive, we can now load the self energy
of the last iteration::

  if previous_present:
    if mpi.is_master_node():
        ar = HDFArchive(dft_filename+'.h5','a')
        S.Sigma_iw << ar['dmft_output']['Sigma_iw']
        del ar
        chemical_potential,dc_imp,dc_energ = SK.load(['chemical_potential','dc_imp','dc_energ'])
    S.Sigma_iw << mpi.bcast(S.Sigma_iw)
    SK.set_mu(chemical_potential)
    SK.set_dc(dc_imp,dc_energ)
    
The self-energy is broadcast from the master node to the slave nodes. Also, the
last saved chemical potential and double counting values are read in and set.

Now we can go to the definition of the self-consistency step. It consists again
of the basic steps discussed in the previous section, with some additional
refinements::

  for iteration_number in range(1,loops+1):
      if mpi.is_master_node(): print "Iteration = ", iteration_number
  
      SK.symm_deg_gf(S.Sigma_iw,orb=0)                        # symmetrise Sigma
      SK.put_Sigma(Sigma_imp = [ S.Sigma_iw ])                # put Sigma into the SumK class
      chemical_potential = SK.calc_mu( precision = prec_mu )  # find the chemical potential for given density
      S.G_iw << SK.extract_G_loc()[0]                         # calc the local Green function
      mpi.report("Total charge of Gloc : %.6f"%S.G_iw.total_density())
  
      # Init the DC term and the real part of Sigma, if no previous runs found:
      if (iteration_number==1 and previous_present==False):
          dm = S.G_iw.density()
          SK.calc_dc(dm, U_interact = U, J_hund = J, orb = 0, use_dc_formula = dc_type)
          S.Sigma_iw << SK.dc_imp[0]['up'][0,0]
  
      # Calculate new G0_iw to input into the solver:
      S.G0_iw << S.Sigma_iw + inverse(S.G_iw)
      S.G0_iw << inverse(S.G0_iw)

      # Solve the impurity problem:
      S.solve(h_int=h_int, **p)
  
      # Solved. Now do post-solution stuff:
      mpi.report("Total charge of impurity problem : %.6f"%S.G_iw.total_density())
  
      # Now mix Sigma and G with factor sigma_mix, if wanted:
      if (iteration_number>1 or previous_present):
          if mpi.is_master_node():
              ar = HDFArchive(dft_filename+'.h5','a')
              mpi.report("Mixing Sigma and G with factor %s"%sigma_mix)
              S.Sigma_iw << sigma_mix * S.Sigma_iw + (1.0-sigma_mix) * ar['dmft_output']['Sigma_iw']
              S.G_iw << sigma_mix * S.G_iw + (1.0-sigma_mix) * ar['dmft_output']['G_iw']
              del ar
          S.G_iw << mpi.bcast(S.G_iw)
          S.Sigma_iw << mpi.bcast(S.Sigma_iw)
  
      # Write the final Sigma and G to the hdf5 archive:
      if mpi.is_master_node():
          ar = HDFArchive(dft_filename+'.h5','a')
          ar['dmft_output']['iterations'] = iteration_number + previous_runs
          ar['dmft_output']['G_0'] = S.G0_iw
          ar['dmft_output']['G_tau'] = S.G_tau
          ar['dmft_output']['G_iw'] = S.G_iw
          ar['dmft_output']['Sigma_iw'] = S.Sigma_iw
          del ar

      # Set the new double counting:
      dm = S.G_iw.density() # compute the density matrix of the impurity problem
      SK.calc_dc(dm, U_interact = U, J_hund = J, orb = 0, use_dc_formula = dc_type)

      # Save stuff into the dft_output group of hdf5 archive in case of rerun:
      SK.save(['chemical_potential','dc_imp','dc_energ'])

This is all we need for the DFT+DMFT calculation. At the end, all results are stored in the hdf5 output file.
