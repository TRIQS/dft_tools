import pytriqs.utility.mpi as mpi
from pytriqs.operators.hamiltonians import *
from pytriqs.archive import HDFArchive
from pytriqs.applications.impurity_solvers.cthyb import *
from pytriqs.gf.local import *
from pytriqs.applications.dft.sumk_dft import *
from pytriqs.applications.dft.converters.wien2k_converter import *

dft_filename='Gd_fcc'
U = 9.6
J = 0.8
beta = 40
loops = 10                       # Number of DMFT sc-loops
sigma_mix = 1.0                  # Mixing factor of Sigma after solution of the AIM
delta_mix = 1.0                  # Mixing factor of Delta as input for the AIM
dc_type = 0                      # DC type: 0 FLL, 1 Held, 2 AMF
use_blocks = True                # use bloc structure from DFT input
prec_mu = 0.0001

# Solver parameters
p = {}
p["max_time"] = -1
p["random_name"] = ""
p["random_seed"] = 123 * mpi.rank + 567
p["length_cycle"] = 50
p["n_warmup_cycles"] = 50
p["n_cycles"] = 5000

Converter = Wien2kConverter(filename=dft_filename, repacking=True)
Converter.convert_dmft_input()
mpi.barrier()

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
# if previous runs are present, no need for recalculating the bloc structure:
calc_blocs = use_blocks and (not previous_present)

SK=SumkDFT(hdf_file=dft_filename+'.h5',use_dft_blocks=calc_blocs)

n_orb = SK.corr_shells[0]['dim']
l = SK.corr_shells[0]['l']
spin_names = ["up","down"]
orb_names = ["%s"%i for i in range(num_orbitals)]
orb_hybridized = False

# Construct U matrix for density-density calculations
gf_struct = set_operator_structure(spin_names,orb_names,orb_hybridized)
# Construct U matrix for density-density calculations
Umat, Upmat = U_matrix_kanamori(n_orb=n_orb, U_int=U, J_hund=J)
# Construct Hamiltonian and solver
H = h_loc_density(spin_names, orb_names, orb_hybridized, U=Umat, Uprime=Upmat, H_dump="H.txt")
S = Solver(beta=beta, gf_struct=gf_struct)

if (previous_present):
  if (mpi.is_master_node()):
      S.Sigma_iw << HDFArchive(dft_filename+'.h5','a')['dmft_output']['Sigma_iw']
  S.Sigma_iw = mpi.bcast(S.Sigma_iw)

for iteration_number in range(1,loops+1):
      if mpi.is_master_node(): print "Iteration = ", i

      SK.symm_deg_gf(S.Sigma_iw,orb=0)                        # symmetrise Sigma
      SK.put_Sigma(Sigma_imp = [ S.Sigma_iw ])                # put Sigma into the SumK class
      chemical_potential = SK.find_mu( precision = prec_mu )  # find the chemical potential for the given density
      S.G_iw << SK.extract_G_loc()[0]                           # extract the local Green function
      mpi.report("Total charge of Gloc : %.6f"%S.G_iw.total_density())

      if ((iteration_number==1)and(previous_present==False)):
          # Init the DC term and the real part of Sigma, if no previous run was found:
          dm = S.G_iw.density()
          SK.set_dc(dm, U_interact = U, J_hund = J, orb = 0, use_dc_formula = dc_type)
          S.Sigma_iw << SK.dc_imp[0]['up'][0,0]

      # now calculate new G0_iw to input into the solver:
      if (mpi.is_master_node()):
          # We can do a mixing of Delta in order to stabilize the DMFT iterations:
          S.G0_iw << S.Sigma_iw + inverse(S.G_iw)
          ar = HDFArchive(dft_filename+'.h5','a')['dmft_output']
          if ((iteration_number>1) or (previous_present)):
              mpi.report("Mixing input Delta with factor %s"%delta_mix)
              Delta = (delta_mix * delta(S.G0_iw)) + (1.0-delta_mix) * ar['Delta_iw']
              S.G0_iw << S.G0_iw + delta(S.G0_iw) - Delta

          ar['Delta_iw'] = delta(S.G0_iw)
          S.G0_iw << inverse(S.G0_iw)
          del ar

      S.G0_iw = mpi.bcast(S.G0_iw)

      # Solve the impurity problem:
      S.solve(h_loc=h_loc, **p)

      # solution done, do the post-processing:
      mpi.report("Total charge of impurity problem : %.6f"%S.G_iw.total_density())

      # Now mix Sigma and G with factor sigma_mix, if wanted:
      if ((iteration_number>1) or (previous_present)):
          if (mpi.is_master_node()):
              ar = HDFArchive(dft_filename+'.h5','a')['dmft_output']
              mpi.report("Mixing Sigma and G with factor %s"%sigma_mix)
              S.Sigma_iw << sigma_mix * S.Sigma_iw + (1.0-sigma_mix) * ar['Sigma_iw']
              S.G_iw << sigma_mix * S.G_iw + (1.0-sigma_mix) * ar['G_iw']
              del ar
          S.G_iw = mpi.bcast(S.G_iw)
          S.Sigma_iw = mpi.bcast(S.Sigma_iw)

      # Write the final Sigma and G to the hdf5 archive:
      if (mpi.is_master_node()):
          ar = HDFArchive(dft_filename+'.h5','a')
          ar['iterations'] = previous_runs + iteration_number
          ar['Sigma_iw'] = S.Sigma_iw
          ar['G_iw'] = S.G_iw
          del ar

      dm = S.G_iw.density() # compute the density matrix of the impurity problem
      # Set the double counting
      SK.set_dc( dm, U_interact = U, J_hund = J, orb = 0, use_dc_formula = dc_type)

      # Save stuff into the dft_output group of hdf5 archive in case of rerun:
      SK.save(['chemical_potential','dc_imp','dc_energ'])

if mpi.is_master_node():
    ar = HDFArchive("dftdmft.h5",'w')
    ar["G_iw"] = S.G_iw
    ar["Sigma_iw"] = S.Sigma_iw
