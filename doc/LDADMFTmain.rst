.. index:: LDA+DMFT calculation

.. _LDADMFTmain:

The LDA+DMFT calculation
========================

After having set up the hdf5 arxive, we can now do our LDA+DMFT calculation. It consists of
initialisation steps, and the actual DMFT self consistency loop.

.. index:: initialisation of LDA+DMFT

Initialisation of the calculation
---------------------------------

Before doing the calculation, we have to intialize all the objects that we will need. The first thing is the 
:class:`SumkLDA` class. It contains all basic routines that are necessary to perform a summation in k-space 
to get the local quantities used in DMFT. It is initialized by::

  from pytriqs.applications.dft.sumk_lda import *
  SK = SumkLDA(hdf_file = filename)

The only necessary parameter is the filename of the hdf5 archive. In addition, there are some optional parameters:

  * `mu`: The chemical potential at initialization. This value is only used, if there is no other value found in the hdf5 arxive. Standard is 0.0
  * `h_field`: External magnetic field, standard is 0.0
  * `use_lda_blocks`: If true, the structure of the density matrix is analysed at initialisation, and non-zero matrix elements 
    are identified. The DMFT calculation is then restricted to these matrix elements, yielding a more efficient solution of the 
    local interaction problem. Also degeneracies in orbital and spin space are recognised, and stored for later use. Standard value is `False`. 
  * `lda_data`, `symm_corr_data`, `par_proj_data`, `symm_par_data`, `bands_data`: These string variables define the subgroups in the hdf5 arxive,
    where the corresponding information is stored. The standard values are consistent with the standard values in :ref:`interfacetowien`.

At initialisation, the necessary data is read from the hdf5 file. If we restart a calculation from a previous one, also the information on
the degenerate shells, the block structure of the density matrix, the chemical potential, and double counting correction are read.

.. index:: Multiband solver

Setting up the Multi-Band Solver
--------------------------------

There is a module that helps setting up the multiband CTQMC solver. It is loaded and initialized by::

  from pytriqs.applications.dft.solver_multiband import *
  S = SolverMultiBand(beta, n_orb, gf_struct = SK.gf_struct_solver[0], map=SK.map[0])

The necessary parameters are the inverse temperature `beta`, the Coulomb interaction `U_interact`, the Hund's rule coupling `J_hund`,
and the number of orbitals `n_orb`. There are again several optional parameters that allow to modify the local Hamiltonian to
specific needs. They are:

  * `gf_struct`: Contains the block structure of the local density matrix. Has to be given in the format as calculated by :class:`SumkLDA`.
  * `map`: If `gf_struct` is given as parameter, also `map` has to be given. This is the mapping from the block structure to a general 
    up/down structure.

The solver method is called later by this statement::

  S.solve(U_interact = U, J_hund = J)

The parameters for the Coulomb interaction `U_interact` and the Hunds coupling `J_hund` are necessary parameters.
The following parameters are optional, by highly recommended to be set:

  * `use_matrix`: If `True`, the interaction matrix is calculated from Slater integrals, which are calculated from `U_interact` and 
    `J_hund`. Otherwise, a Kanamori representation is used. Attention: We define the intraorbital interaction as 
    `U_interact+2J_hund`, the interorbital interaction for opposite spins as `U_interact`, and interorbital for equal spins as 
    `U_interact-J_hund`!
  * `T`: A matrix that transforms the interaction matrix from spherical harmonics, to a symmetry adapted basis. Only effective, if 
    `use_matrix=True`.
  * `l`: Orbital quantum number. Again, only effective for Slater parametrisation.
  * `deg_shells`: A list that gives the degeneracies of the orbitals. It is used to set up a global move of the CTQMC solver.
  * `use_spinflip`: If `True`, the full rotationally-invariant interaction is used. Otherwise, only density-density terms are
    kept in the local Hamiltonian.
  * `dim_reps`: If only a subset of the full d-shell is used a correlated orbtials, one can specify here the dimensions of all the subspaces
    of the d-shell, i.e. t2g and eg. Only effective for Slater parametrisation.
  * `irep`: The index in the list `dim_reps` of the subset that is used. Only effective for Slater parametrisation.

Most of above parameters can be taken directly from the :class:`SumkLDA` class, without defining them by hand. We will see a specific example 
at the end of this tutorial.

After initialisation, several other CTQMC parameters can be set (see CTQMC doc). The most important are:

  * `S.n_cycles`: Number of QMC cycles per node.
  * `S.n_warmup_cycles`: Number of iterations used for thermalisation




.. index:: LDA+DMFT loop, one-shot calculation

Doing the DMFT loop
-------------------

Having initialised the SumK class and the Solver, we can proceed with the DMFT loop itself. As explained in the tutorial, we have to 
set up the loop over DMFT iterations and the self-consistency condition::

  n_loops = 5
  for iteration_number in range(n_loops) :            # start the DMFT loop

          SK.put_Sigma(Sigma_imp = [ S.Sigma ])      # Put self energy to the SumK class
          chemical_potential = SK.find_mu()          # find the chemical potential for the given density
          S.G <<= SK.extract_G_loc()[0]              # extract the local Green function
          S.G0 <<= inverse(S.Sigma + inverse(S.G))   # finally get G0, the input for the Solver

          S.Solve(U_interact = U, J_hund = J)                                  # now solve the impurity problem

	  dm = S.G.density()                         # density matrix of the impurity problem  
          SK.set_dc( dm, U_interact = U, J_hund = J, use_dc_formula = 0)     # Set the double counting term
          SK.save()                                  # save everything to the hdf5 arxive

These basic steps are enough to set up the basic DMFT Loop. For a detailed description of the :class:`SumkLDA` routines,
see the reference manual. After the self-consistency steps, the solution of the Anderson impurity problem is calculation by CTQMC. 
Different to model calculations, we have to do a few more steps after this, because of the double-counting correction. We first 
calculate the density of the impurity problem. Then, the routine `set_dc` takes as parameters this density matrix, the 
Coulomb interaction, Hund's rule coupling, and the type of double-counting that should be used. Possible values for `use_dc_formula` are:

  * `0`: Full-localised limit
  * `1`: DC formula as given in K. Held, Adv. Phys. 56, 829 (2007).
  * `2`: Around-mean-field

At the end of the calculation, we can save the Greens function and self energy into a file::

  from pytriqs.archive import HDFArchive
  import pytriqs.utility.mpi as mpi
  if mpi.is_master_node():
      R = HDFArchive("YourLDADMFTcalculation.h5",'w')
      R["G"] = S.G
      R["Sigma"] = S.Sigma

This is it! 

These are the essential steps to do a one-shot LDA+DMFT calculation. For full charge-self consistent calculations, there are some more things
to consider, which we will see later on.
