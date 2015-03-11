.. index:: DFT+DMFT calculation

.. highlight:: python

.. _DFTDMFTmain:

The DFT+DMFT calculation
========================

After having set up the hdf5 archive, we can now do our DFT+DMFT calculation. It consists of
initialisation steps, and the actual DMFT self consistency loop.

.. index:: initialisation of DFT+DMFT

Initialisation of the calculation
---------------------------------

Before doing the calculation, we have to intialize all the objects that we will need. The first thing is the 
:class:`SumkDFT` class. It contains all basic routines that are necessary to perform a summation in k-space 
to get the local quantities used in DMFT. It is initialized by::

  from pytriqs.applications.dft.sumk_dft import *
  SK = SumkDFT(hdf_file = filename)


This is the reference for the function:

.. function::  SumkDFT(hdf_file, h_field = 0.0, use_dft_blocks = False, dft_data = 'dft_input', symmcorr_data = 'dft_symmcorr_input', parproj_data = 'dft_parproj_input', symmpar_data = 'dft_symmpar_input', bands_data = 'dft_bands_input', transp_data = 'dft_transp_input', misc_data = 'dft_misc_input')

The parameters needed to initialise SumkDFT as as follows (with any default variables as shown above):

=========================   ====================  ===========================================================================
Name                        Type                  Meaning
=========================   ====================  ===========================================================================
hdf_file                    String                Name of the main hdf5 file containing converted dft information. 
h_field                     Scalar                External magnetic field.
use_dft_blocks              Boolean               Analyse the structure of the density matrix at initialisation,
                                                  and identify non-zero matrix elements. 
                                                  The DMFT calculation is then restricted to these matrix elements, 
                                                  yielding a more efficient solution of the local interaction problem. 
                                                  Degeneracies in orbital and spin space are also identified and stored for later use.
dft_data                    String                hdf5 subgroup containing required DFT data.
symmcorr_data               String                hdf5 subgroup containing correlated shell symmetry data.                    
parproj_data                String                hdf5 subgroup containing partial projector data.
symmpar_data                String                hdf5 subgroup containing non-correlated shell symmetry data. 
bands_data                  String                hdf5 subgroup containing DFT band data.
transp_data                 String                hdf5 subgroup containing optics/transport data.
misc_data                   String                hdf5 subgroup containing miscellaneous DFT data.
=========================   ====================  ===========================================================================

.. index:: Multiband solver

Setting up the Multi-Band Solver
--------------------------------

The next step is to setup the impurity solver. 
For more details here, see the `CTHYB <http://ipht.cea.fr/triqs/1.2/applications/cthyb/>`_ documentation.

.. This is initialised as follows::
.. 
..     from pytriqs.applications.impurity_solvers.cthyb import *
..     beta = 40.0
..     gf_struct = SK.gf_struct_solver[0]
..     S = Solver(beta=beta, gf_struct=gf_struct)
.. 
.. The solver method is called later by this statement::
.. 
..     S.solve(h_loc=h_loc, **p)
.. 
.. where `p` represents the solve parameters.
.. 
.. There is a module that helps setting up the multiband CTQMC solver. It is loaded and initialized by::
.. 
..   from pytriqs.applications.dft.solver_multiband import *
..   S = SolverMultiBand(beta, n_orb, gf_struct = SK.gf_struct_solver[0], map=SK.map[0])
.. 
.. The necessary parameters are the inverse temperature `beta`, the Coulomb interaction `U_interact`, the Hund's rule coupling `J_hund`,
.. and the number of orbitals `n_orb`. There are again several optional parameters that allow the tailoring of the local Hamiltonian to
.. specific needs. They are:
.. 
..   * `gf_struct`: The block structure of the local density matrix given in the format calculated by :class:`SumkDFT`.
..   * `map`: If `gf_struct` is given as parameter, `map` also must be given. This is the mapping from the block structure to a general 
..     up/down structure.
.. 
.. The solver method is called later by this statement::
.. 
..   S.solve(U_interact,J_hund,use_spinflip=False,use_matrix=True,
..                    l=2,T=None, dim_reps=None, irep=None, n_cycles =10000,
..                    length_cycle=200,n_warmup_cycles=1000)
.. 
.. The parameters for the Coulomb interaction `U_interact` and the Hund's coupling `J_hund` are necessary input parameters. The rest are optional 
.. parameters for which default values are set. Generally, they should be reset for the problem at hand. Here is a description of the parameters:
.. 
..   * `use_matrix`: If `True`, the interaction matrix is calculated from Slater integrals, which are computed from `U_interact` and 
..     `J_hund`. Otherwise, a Kanamori representation is used. Attention: We define the intraorbital interaction as 
..     `U_interact`, the interorbital interaction for opposite spins as `U_interact-2*J_hund`, and interorbital for equal spins as 
..     `U_interact-3*J_hund`.
..   * `T`: The matrix that transforms the interaction matrix from spherical harmonics to a symmetry-adapted basis. Only effective for Slater
..      parametrisation, i.e. `use_matrix=True`.
..   * `l`: The orbital quantum number. Again, only effective for Slater parametrisation, i.e. `use_matrix=True`.
..   * `use_spinflip`: If `True`, the full rotationally-invariant interaction is used. Otherwise, only density-density terms are
..     kept in the local Hamiltonian.
..   * `dim_reps`: If only a subset of the full d-shell is used as correlated orbtials, one can specify here the dimensions of all the subspaces
..     of the d-shell, i.e. t2g and eg. Only effective for Slater parametrisation.
..   * `irep`: The index in the list `dim_reps` of the subset that is used. Only effective for Slater parametrisation.
..   * `n_cycles`: Number of CTQMC cycles (a sequence of moves followed by a measurement) per core. The default value of 10000 is the minimum, and generally should be increased.
..   * `length_cycle`: Number of CTQMC moves per one cycle.
..   * `n_warmup_cycles`: Number of initial CTQMC cycles before measurements start. Usually of order of 10000, sometimes needs to be increased significantly.
.. 
.. Most of above parameters can be taken directly from the :class:`SumkDFT` class, without defining them by hand. We will see a specific example 
.. at the end of this tutorial.


.. index:: DFT+DMFT loop, one-shot calculation

Doing the DMFT loop
-------------------

Having initialised the SumK class and the Solver, we can proceed with the DMFT loop itself. As explained in the tutorial, we have to 
set up the loop over DMFT iterations and the self-consistency condition::

  n_loops = 5
  for iteration_number in range(n_loops) :            # start the DMFT loop

          SK.put_Sigma(Sigma_imp = [ S.Sigma ])              # Put self energy to the SumK class
          chemical_potential = SK.calc_mu()                  # calculate the chemical potential for the given density
          S.G_iw << SK.extract_G_loc()[0]                    # extract the local Green function
          S.G0_iw << inverse(S.Sigma_iw + inverse(S.G_iw))   # finally get G0, the input for the Solver

          S.solve(h_loc=h_loc, **p)                          # now solve the impurity problem

	  dm = S.G_iw.density()                                                 # Density matrix of the impurity problem  
          SK.calc_dc(dm, U_interact=U, J_hund=J, orb=0, use_dc_formula=dc_type) # Set the double counting term
          SK.save(['chemical_potential','dc_imp','dc_energ'])                   # Save data in the hdf5 archive

These basic steps are enough to set up the basic DMFT Loop. For a detailed description of the :class:`SumkDFT` routines,
see the reference manual. After the self-consistency steps, the solution of the Anderson impurity problem is calculation by CTQMC. 
Different to model calculations, we have to do a few more steps after this, because of the double-counting correction. We first 
calculate the density of the impurity problem. Then, the routine `calc_dc` takes as parameters this density matrix, the 
Coulomb interaction, Hund's rule coupling, and the type of double-counting that should be used. Possible values for `use_dc_formula` are:

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

These are the essential steps to do a one-shot DFT+DMFT calculation. For full charge-self consistent calculations, there are some more things
to consider, which we will see later on.
