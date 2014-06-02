.. _advanced:

A more advanced example
=======================

Normally, one wants to adjust some more parameters in order to make the calculation more efficient. Here, we
will see a more advanced example, which is also suited for parallel execution. 
First, we load the necessary modules::

  from pytriqs.applications.dft.sumk_lda import *
  from pytriqs.applications.dft.converters.wien2k_converter import *
  from pytriqs.applications.dft.solver_multiband import *
  from pytriqs.gf.local import *
  from pytriqs.archive import *

Then we define some parameters::

  lda_filename='srvo3'
  U = 2.7
  J = 0.65
  beta = 40
  loops =  10                     # Number of DMFT sc-loops
  mix = 0.8                        # Mixing factor of Sigma after solution of the AIM
  Delta_mix = 1.0                  # Mixing factor of Delta as input for the AIM
  dc_type = 1                      # DC type: 0 FLL, 1 Held, 2 AMF
  use_blocks = True                # use bloc structure from LDA input
  use_matrix = False               # True: Slater parameters, False: Kanamori parameters U+2J, U, U-J
  use_spinflip = False             # use the full rotational invariant interaction?
  prec_mu = 0.0001
  qmc_cycles = 20000
  length_cycle = 200
  warming_iterations = 2000


Most of these parameters are self-explaining. The first, `lda_filename`, gives the filename of the input files. 
The next step, as described in the previous section, is to convert the input files::

  Converter = Wien2kConverter(filename=lda_filename, repacking=True)
  Converter.convert_dmft_input()
  mpi.barrier()

The command ``mpi.barrier()`` ensures that all nodes wait until the conversion of the input is finished on the master
node. After the conversion, we can check in the hdf5 archive, if previous runs are present, or if we have to start
from scratch::

  previous_runs = 0
  previous_present = False
  if mpi.is_master_node():
      ar = HDFArchive(lda_filename+'.h5','a')
      if 'iterations' in ar:
          previous_present = True
          previous_runs = ar['iterations']
      del ar
  previous_runs    = mpi.bcast(previous_runs)
  previous_present = mpi.bcast(previous_present)
  # if previous runs are present, no need for recalculating the bloc structure:
  calc_blocs = use_blocks and (not previous_present)

Now we can use all this information to initialise the :class:`SumkLDA` class::

  SK=SumkLDA(hdf_file=lda_filename+'.h5',use_lda_blocks=calc_blocs)

If there was a previous run, we know already about the block structure, and therefore `UseLDABlocs` is set to `False`.
The next step is to initialise the Solver::

  Norb = SK.corr_shells[0][3]
  l = SK.corr_shells[0][2]
  S = SolverMultiBand(beta=beta,n_orb=Norb,gf_struct=SK.gf_struct_solver[0],map=SK.map[0])

As we can see, many options of the solver are set by properties of the :class:`SumkLDA` class, so we don't have
to set them manually. 

If there are previous runs stored in the hdf5 archive, we can now load the self energy
of the last iteration::

  if (previous_present):
    if (mpi.is_master_node()):
        ar = HDFArchive(lda_filename+'.h5','a')
        S.Sigma <<= ar['SigmaF']
        del ar
    S.Sigma = mpi.bcast(S.Sigma)
    
The last command is the broadcasting of the self energy from the master node to the slave nodes. 
Now we can go to the definition of the self-consistency step. It consists again of the basic steps discussed in the 
previous section, with some additional refinement::

  for iteration_number in range(1,loops+1) :
  
        SK.symm_deg_gf(S.Sigma,orb=0)                           # symmetrise Sigma
        SK.put_Sigma(Sigma_imp = [ S.Sigma ])                   # put Sigma into the SumK class:
  
        chemical_potential = SK.find_mu( precision = prec_mu )  # find the chemical potential
        S.G <<= SK.extract_G_loc()[0]                           # calculation of the local Green function
        mpi.report("Total charge of Gloc : %.6f"%S.G.total_density())
  
        if ((iteration_number==1)and(previous_present==False)):
            # Init the DC term and the real part of Sigma, if no previous run was found:
            dm = S.G.density()
            SK.set_dc( dm, U_interact = U, J_hund = J, orb = 0, use_dc_formula = dc_type)
            S.Sigma <<= SK.dc_imp[0]['up'][0,0]
  
        S.G0 <<= inverse(S.Sigma + inverse(S.G))
  
        # Solve the impurity problem:
        S.solve(U_interact=U,J_hund=J,use_spinflip=use_spinflip,use_matrix=use_matrix,
                     l=l,T=SK.T[0], dim_reps=SK.dim_reps[0], irep=2, deg_orbs=SK.deg_shells[0],n_cycles =qmc_cycles,
                     length_cycle=length_cycle,n_warmup_cycles=warming_iterations)
  
        # solution done, do the post-processing:
        mpi.report("Total charge of impurity problem : %.6f"%S.G.total_density())
  
        S.Sigma <<=(inverse(S.G0)-inverse(S.G))
        # Solve the impurity problem:
        S.solve(U_interact=U,J_hund=J,use_spinflip=use_spinflip,use_matrix=use_matrix,
                     l=l,T=SK.T[0], dim_reps=SK.dim_reps[0], irep=2, deg_orbs=SK.deg_shells[0],n_cycles =qmc_cycles,
                     length_cycle=length_cycle,n_warmup_cycles=warming_iterations)
  
        # solution done, do the post-processing:
        mpi.report("Total charge of impurity problem : %.6f"%S.G.total_density())
  
        S.Sigma <<=(inverse(S.G0)-inverse(S.G))
  
        # Now mix Sigma and G with factor Mix, if wanted:
        if ((iteration_number>1) or (previous_present)):
            if (mpi.is_master_node()):
                ar = HDFArchive(lda_filename+'.h5','a')
                mpi.report("Mixing Sigma and G with factor %s"%mix)
                S.Sigma <<= mix * S.Sigma + (1.0-mix) * ar['Sigma']
                S.G <<= mix * S.G + (1.0-mix) * ar['GF']
                del ar
            S.G = mpi.bcast(S.G)
            S.Sigma = mpi.bcast(S.Sigma)
  
        # Write the final Sigma and G to the hdf5 archive:
        if (mpi.is_master_node()):
            ar = HDFArchive(lda_filename+'.h5','a')
            ar['iterations'] = previous_runs + iteration_number
            ar['Sigma'] = S.Sigma
            ar['GF'] = S.G
            del ar

        # Now set new double counting:
        dm = S.G.density()
        SK.set_dc( dm, U_interact = U, J_hund = J, orb = 0, use_dc_formula = dc_type)
  
        #Save stuff:
        SK.save()
                                
This is all we need for the LDA+DMFT calculation. At the end, all results are stored in the hdf5 output file.



