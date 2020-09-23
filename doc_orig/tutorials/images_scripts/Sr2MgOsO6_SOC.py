# Import the modules:
from triqs_dft_tools.sumk_dft import *
from triqs.gf import *
from h5 import HDFArchive
from triqs.operators.util import *
from triqs.operators.util.U_matrix import *
from triqs_cthyb import *
import triqs.utility.mpi as mpi

# Init the SumK class:
filename = 'Sr2MgOsO6_SOC.h5'
SK = SumkDFT(hdf_file=filename,use_dft_blocks=True)

# Find diagonal local basis set:
SK.calculate_diagonalization_matrix(prop_to_be_diagonal='eal',calc_in_solver_blocks=True)

###########################
# Now we pick the orbitals:
# BE CAREFUL: THIS NEEDS TO BE DONE PROPERLY 
# AND IS DIFFERENT FORM CASE TO CASE!
SK.block_structure.pick_gf_struct_solver([{'ud_0': [0,1,2],'ud_1': [0,1,2]}])
###########################

# Now we set up the U matrix, first in cubic (Wien2k) convention:
U = 2.0
J = 0.2
U_sph = U_matrix(l=2, U_int=U, J_hund=J)
U_sph = np.kron(np.reshape(np.eye(2),(1,2,1,2)),np.kron(np.reshape(np.eye(2),(2,1,2,1)),U_sph))
U_mat = transform_U_matrix(U_sph, SK.T[0].conjugate())

# Now we set up the interaction Hamiltonian
h_sumk = h_int_slater(['ud'], range(10), U_mat,  off_diag=True, complex=True)
# And convert it to the solver structure
h_int = SK.block_structure.convert_operator(h_sumk)

# Solver Init:
beta = 40.0
S = Solver(beta=beta, gf_struct=SK.block_structure.gf_struct_solver_list[0])

# Solver parameters:
p = {}
# solver
p["random_seed"] = 123 * mpi.rank + 567
p["length_cycle"] = 200
p["n_warmup_cycles"] = 100000
p["n_cycles"] = 3000000
# tail fit
p["perform_tail_fit"] = True
p["fit_max_moment"] = 4
p["fit_min_w"] = 4.0
p["fit_max_w"] = 8.0

# double counting correction:
dc_type = 0  # FLL
# DMFT loops:
n_loops = 1

#for first iteration, add the outout group:
if mpi.is_master_node():
    with HDFArchive(filename) as ar:
        if (not ar.is_group('dmft_output')):
            ar.create_group('dmft_output')

for iteration_number in range(1,n_loops+1):

    mpi.report("Iteration = %s"%iteration_number)
    
    SK.set_Sigma([ S.Sigma_iw ])                    # put Sigma into the SumK class
    chemical_potential = SK.calc_mu( precision = 0.01 )  # find the chemical potential for given density
    S.G_iw << SK.extract_G_loc()[0]

    if (iteration_number==1):
        # Put Hartree energy on Re Sigma
        dm = S.G_iw.density()
        SK.calc_dc(dm, U_interact = U, J_hund = J, orb = 0, use_dc_formula = dc_type)
        S.Sigma_iw << SK.block_structure.convert_matrix(SK.dc_imp[0],space_from='sumk',space_to='solver')['ud_0'][0,0]

    mpi.report("Orbital densities of local Green function :")
    for s,gf in S.G_iw:
        mpi.report("Orbital %s: %s"%(s,dm[s].real))
    mpi.report("Total charge of Gloc : %.6f"%S.G_iw.total_density().real)

    # Calculate new G0_iw to input into the solver:
    S.G0_iw << S.Sigma_iw + inverse(S.G_iw)
    S.G0_iw << inverse(S.G0_iw)

    # Solve the impurity problem:
    S.solve(h_int=h_int, **p)

    # Solved. Now do post-solution stuff:
    dm = S.G_iw.density()
    mpi.report("Orbital densities of impurity Green function:")
    for s,gf in S.G_iw:
        mpi.report("Orbital %s: %s"%(s,dm[s].real))
    mpi.report("Total charge of impurity problem : %.6f"%S.G_iw.total_density().real)

    # Write the final Sigma and G to the hdf5 archive:
    if mpi.is_master_node():
        with HDFArchive(filename) as ar:
              ar['dmft_output']['iterations'] = iteration_number
              ar['dmft_output']['G_0'] = S.G0_iw
              ar['dmft_output']['G_tau'] = S.G_tau
              ar['dmft_output']['G_iw'] = S.G_iw
              ar['dmft_output']['Sigma_iw_%s'%iteration_number] = S.Sigma_iw

    # Set the new double counting:
    dm = S.G_iw.density() # compute the density matrix of the impurity problem
    SK.calc_dc(dm, U_interact = U, J_hund = J, orb = 0, use_dc_formula = dc_type)

    # Save stuff into the user_data group of hdf5 archive in case of rerun:
    SK.save(['chemical_potential','dc_imp','dc_energ'])
    
