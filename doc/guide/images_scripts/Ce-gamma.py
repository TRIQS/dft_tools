from pytriqs.applications.dft.sumk_dft import *
from pytriqs.applications.dft.converters.wien2k_converter import *
from pytriqs.applications.impurity_solvers.hubbard_I.hubbard_solver import Solver


import os
dft_filename = os.getcwd().rpartition('/')[2]

beta = 40
U_int = 6.00
J_hund = 0.70
Loops =  5                       # Number of DMFT sc-loops
Mix = 0.7                        # Mixing factor in QMC
DC_type = 0                      # 0...FLL, 1...Held, 2... AMF, 3...Lichtenstein
chemical_potential_init=0.0      # initial chemical potential

HDFfilename = dft_filename+'.h5'

# Convert DMFT input:
Converter = Wien2kConverter(filename=filename)
Converter.convert_dft_input()
mpi.barrier()

#check if there are previous runs:
previous_runs = 0
previous_present = False
if mpi.is_master_node():
    f = HDFArchive(filename+'.h5','a')
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

# Init the SumK class
SK=SumkDFT(hdf_file=dft_filename+'.h5',use_dft_blocks=False)

Norb = SK.corr_shells[0]['dim']
l    = SK.corr_shells[0]['l']

# Init the Hubbard-I solver:
S = Solver(beta = beta, l = l)

chemical_potential=chemical_potential_init
# load previous data: old self-energy, chemical potential, DC correction
if previous_present:
    if mpi.is_master_node():
        ar = HDFArchive(dft_filename+'.h5','a')
        S.Sigma << ar['dmft_output']['Sigma']
        del ar
        chemical_potential,dc_imp,dc_energ = SK.load(['chemical_potential','dc_imp','dc_energ'])
    S.Sigma << mpi.bcast(S.Sigma)
    SK.set_mu(chemical_potential)
    SK.set_dc(dc_imp,dc_energ)
    

# DMFT loop:
for iteration_number in range(1,Loops+1):
    
        itn = iteration_number + previous_runs
       
        # put Sigma into the SumK class:
        SK.put_Sigma(Sigma_imp = [ S.Sigma ])

        # Compute the SumK, possibly fixing mu by dichotomy
        chemical_potential = SK.calc_mu( precision = 0.000001 )
                  
        # Density:
        S.G <<= SK.extract_G_loc()[0]
        mpi.report("Total charge of Gloc : %.6f"%S.G.total_density())

        # calculated DC at the first run to have reasonable initial non-interacting atomic level positions
        if ((iteration_number==1)and(previous_present==False)):
            dc_value_init=U_int/2.0
            dm=S.G.density()
	        SK.calc_dc( dm, U_interact = U_int, J_hund = J_hund, orb = 0, use_dc_formula = DC_type, use_dc_value=dc_value_init)

        # calculate non-interacting atomic level positions:
        eal = SK.eff_atomic_levels()[0]
        S.set_atomic_levels( eal = eal )

        # solve it:
        S.solve(U_int = U_int, J_hund = J_hund, verbosity = 1)

        # Now mix Sigma and G with factor Mix, if wanted:
        if (iteration_number>1 or previous_present):
            if (mpi.is_master_node() and (sigma_mix<1.0)):
                ar = HDFArchive(dft_filename+'.h5','a')
                mpi.report("Mixing Sigma and G with factor %s"%sigma_mix)
                S.Sigma << sigma_mix * S.Sigma + (1.0-sigma_mix) * ar['dmft_output']['Sigma']
                S.G << sigma_mix * S.G + (1.0-sigma_mix) * ar['dmft_output']['G']
                del ar
            S.G << mpi.bcast(S.G)
            S.Sigma << mpi.bcast(S.Sigma)
     
        
        # after the Solver has finished, set new double counting: 
        dm = S.G.density()
        SK.calc_dc( dm, U_interact = U_int, J_hund = J_hund, orb = 0, use_dc_formula = DC_type )

        # correlation energy calculations:
        correnerg = 0.5 * (S.G * S.Sigma).total_density()
        mpi.report("Corr. energy = %s"%correnerg)

        # store the impurity self-energy, GF as well as correlation energy in h5
        if mpi.is_master_node():
            ar = HDFArchive(dft_filename+'.h5','a')
            ar['dmft_output']['iterations'] = iteration_number + previous_runs
            ar['dmft_output']['G'] = S.G
            ar['dmft_output']['Sigma'] = S.Sigma
            del ar

        #Save essential SumkDFT data:
        SK.save(['chemical_potential','dc_imp','dc_energ','correnerg'])
        if (mpi.is_master_node()):
            print 'DC after solver: ',SK.dc_imp[0]

        # print out occupancy matrix of Ce 4f
        mpi.report("Orbital densities of impurity Green function:")
        for s in dm:
            mpi.report("Block %s: "%s)
            for ii in range(len(dm[s])):
                str = ''
                for jj in range(len(dm[s])):
                    if (dm[s][ii,jj].real>0):
                        str += "   %.4f"%(dm[s][ii,jj].real)
                    else:
                        str += "  %.4f"%(dm[s][ii,jj].real)
                mpi.report(str)
        mpi.report("Total charge of impurity problem : %.6f"%S.G.total_density())


# find exact chemical potential
SK.chemical_potential = SK.calc_mu( precision = 0.000001 )

# calculate and save occupancy matrix in the Bloch basis for Wien2k charge denity recalculation
dN,d = SK.calc_density_correction(filename = dft_filename+'.qdmft')

mpi.report("Trace of Density Matrix: %s"%d)

# store correlation energy contribution to be read by Wien2ki and then included to DFT+DMFT total energy
if (mpi.is_master_node()):
    correnerg -= DCenerg[0]
    f=open(dft_filename+'.qdmft','a')
    f.write("%.16f\n"%correnerg)
    f.close()
