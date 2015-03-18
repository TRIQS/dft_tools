from pytriqs.applications.dft.sumk_dft import *
from pytriqs.applications.dft.converters.wien2k_converter import *
from pytriqs.applications.impurity_solvers.hubbard_I.hubbard_solver import Solver

lda_filename = 'Ce-gamma'
beta = 40
U_int = 6.00
J_hund = 0.70
Loops =  5                       # Number of DMFT sc-loops
Mix = 0.7                        # Mixing factor in QMC
DC_type = 0                      # 0...FLL, 1...Held, 2... AMF, 3...Lichtenstein
DC_Mix = 1.0                     # 1.0 ... all from imp; 0.0 ... all from Gloc   
useBlocs = False                 # use bloc structure from LDA input
useMatrix = True                 # use the U matrix calculated from Slater coefficients instead of (U+2J, U, U-J)
chemical_potential_init=0.0      # initial chemical potential

HDFfilename = lda_filename+'.h5'

# Convert DMFT input:
# Can be commented after the first run
Converter = Wien2kConverter(filename=lda_filename)
Converter.convert_dft_input()

#check if there are previous runs:
previous_runs = 0
previous_present = False

if mpi.is_master_node():
    ar = HDFArchive(HDFfilename,'a')
    if 'iterations' in ar:
        previous_present = True
        previous_runs = ar['iterations']
    else: 
        previous_runs = 0
        previous_present = False
    del ar

mpi.barrier()
previous_runs    = mpi.bcast(previous_runs)
previous_present = mpi.bcast(previous_present)

# Init the SumK class
SK=SumkDFT(hdf_file=lda_filename+'.h5',use_dft_blocks=False)

Norb = SK.corr_shells[0]['dim']
l    = SK.corr_shells[0]['l']

# Init the Hubbard-I solver:
S = Solver(beta = beta, l = l)

chemical_potential=chemical_potential_init
# load previous data: old self-energy, chemical potential, DC correction
if (previous_present):
    mpi.report("Using stored data for initialisation")
    if (mpi.is_master_node()):
        ar = HDFArchive(HDFfilename,'a')
        S.Sigma <<= ar['SigmaF']
        del ar
        things_to_load=['chemical_potential','dc_imp']
        old_data=SK.load(things_to_load)
        chemical_potential=old_data[0]
        SK.dc_imp=old_data[1]
    S.Sigma = mpi.bcast(S.Sigma)
    chemical_potential=mpi.bcast(chemical_potential)
    SK.dc_imp=mpi.bcast(SK.dc_imp)


# DMFT loop:
for Iteration_Number in range(1,Loops+1):
    
        itn = Iteration_Number + previous_runs
       
        # put Sigma into the SumK class:
        SK.put_Sigma(Sigma_imp = [ S.Sigma ])

        # Compute the SumK, possibly fixing mu by dichotomy
        if SK.density_required and (Iteration_Number > 1):
            chemical_potential = SK.calc_mu( precision = 0.000001 )
        else:
            mpi.report("No adjustment of chemical potential\nTotal density  = %.3f"%SK.total_density(mu=chemical_potential))

        # Density:
        S.G <<= SK.extract_G_loc()[0]
        mpi.report("Total charge of Gloc : %.6f"%S.G.total_density())

        # calculated DC at the first run to have reasonable initial non-interacting atomic level positions
        if ((Iteration_Number==1)and(previous_present==False)):
            dc_value_init=U_int/2.0
            dm=S.G.density()
	    SK.calc_dc( dm, U_interact = U_int, J_hund = J_hund, orb = 0, use_dc_formula = DC_type, use_dc_value=dc_value_init)

        # calculate non-interacting atomic level positions:
        eal = SK.eff_atomic_levels()[0]
        S.set_atomic_levels( eal = eal )

        # solve it:
        S.solve(U_int = U_int, J_hund = J_hund, verbosity = 1)

        # Now mix Sigma and G:
        if ((itn>1)or(previous_present)):
            if (mpi.is_master_node()and (Mix<1.0)):
                ar = HDFArchive(HDFfilename,'r')
                mpi.report("Mixing Sigma and G with factor %s"%Mix)
                if ('SigmaF' in ar):
                    S.Sigma <<= Mix * S.Sigma + (1.0-Mix) * ar['SigmaF']
                if ('GF' in ar):
                    S.G <<= Mix * S.G + (1.0-Mix) * ar['GF']
                del ar
            S.G = mpi.bcast(S.G)
            S.Sigma = mpi.bcast(S.Sigma)
        
        # after the Solver has finished, set new double counting: 
        dm = S.G.density()
        SK.calc_dc( dm, U_interact = U_int, J_hund = J_hund, orb = 0, use_dc_formula = DC_type )

        # correlation energy calculations:
        correnerg = 0.5 * (S.G * S.Sigma).total_density()
        mpi.report("Corr. energy = %s"%correnerg)

        # store the impurity self-energy, GF as well as correlation energy in h5
        if (mpi.is_master_node()):
            ar = HDFArchive(HDFfilename,'a')
            ar['iterations'] = itn
            ar['chemical_cotential%s'%itn] = chemical_potential
            ar['SigmaF'] = S.Sigma
            ar['GF'] = S.G
            ar['correnerg%s'%itn] = correnerg
            ar['DCenerg%s'%itn] = SK.dc_energ
            del ar

        #Save essential SumkDFT data:
        things_to_save=['chemical_potential','dc_energ','dc_imp']
        SK.save(things_to_save)
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
if (SK.density_required):
    SK.chemical_potential = SK.calc_mu( precision = 0.000001 )

# calculate and save occupancy matrix in the Bloch basis for Wien2k charge denity recalculation
dN,d = SK.calc_density_correction(filename = lda_filename+'.qdmft')

mpi.report("Trace of Density Matrix: %s"%d)

# store correlation energy contribution to be read by Wien2ki and then included to DFT+DMFT total energy
if (mpi.is_master_node()):
    ar = HDFArchive(HDFfilename)
    itn = ar['iterations'] 
    correnerg = ar['correnerg%s'%itn] 
    DCenerg = ar['DCenerg%s'%itn]
    del ar
    correnerg -= DCenerg[0]
    f=open(lda_filename+'.qdmft','a')
    f.write("%.16f\n"%correnerg)
    f.close()
