from pytriqs.applications.dft.sumk_dft_tools import *
from pytriqs.applications.dft.converters.wien2k_converter import *
from pytriqs.applications.impurity_solvers.hubbard_I.hubbard_solver import Solver

# Creates the data directory, cd into it:
#Prepare_Run_Directory(DirectoryName = "Ce-Gamma") 
dft_filename = 'Ce-gamma'
Beta =  40
U_int = 6.00
J_hund = 0.70
DC_type = 0                      # 0...FLL, 1...Held, 2... AMF, 3...Lichtenstein
load_previous = True              # load previous results
useBlocs = False                 # use bloc structure from DFT input
useMatrix = True                 # use the U matrix calculated from Slater coefficients instead of (U+2J, U, U-J)
ommin=-4.0
ommax=6.0
N_om=2001
broadening = 0.02

HDFfilename = dft_filename+'.h5'

# Convert DMFT input:
# Can be commented after the first run
Converter = Wien2kConverter(filename=dft_filename,repacking=True)
Converter.convert_dmft_input()
Converter.convert_parproj_input()

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

# if previous runs are present, no need for recalculating the bloc structure
# It has to be commented, if you run this script for the first time, starting
# from a converted h5 archive.

# Init the SumK class
SK = SumkDFTTools(hdf_file=dft_filename+'.h5',use_dft_blocks=False)


if (mpi.is_master_node()):
    print 'DC after reading SK: ',SK.dc_imp[SK.invshellmap[0]]

N = SK.corr_shells[0]['dim']
l = SK.corr_shells[0]['l']

# Init the Solver:
S = Solver(beta = Beta, l = l)
S.Nmoments= 8

# set atomic levels:
eal = SK.eff_atomic_levels()[0]
S.set_atomic_levels( eal = eal )
S.GF_realomega(ommin=ommin, ommax = ommax, N_om=N_om,U_int=U_int,J_hund=J_hund)
SK.put_Sigma(Sigma_imp = [S.Sigma])
SK.dos_partial(broadening=broadening)
