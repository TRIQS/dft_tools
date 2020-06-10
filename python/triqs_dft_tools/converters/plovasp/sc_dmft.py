 
################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2011 by M. Ferrero, O. Parcollet
#
# DFT tools: Copyright (C) 2011 by M. Aichhorn, L. Pourovskii, V. Vildosola
#
# PLOVasp: Copyright (C) 2015 by O. E. Peil
#
# TRIQS is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# TRIQS. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

import errno
import os
import re
import time
import signal
import sys
import triqs.utility.mpi as mpi
from . import converter
from shutil import copyfile

xch = sys.excepthook
def excepthook(typ, value, traceback):
    xch(typ, value, traceback)
    if mpi.MPI.COMM_WORLD.size > 1:
        mpi.MPI.COMM_WORLD.Abort(1)
sys.excepthook = excepthook

debug = False
#
# Helper functions
#
def sigint_handler(signal, frame):
    raise SystemExit(1)

def is_vasp_lock_present():
    res_bool = False
    if mpi.is_master_node():
        res_bool = os.path.isfile('./vasp.lock')
    res_bool = mpi.bcast(res_bool)
    return res_bool

def is_vasp_running(vasp_pid):
    """
    Tests if VASP initial process is still alive.
    """
    pid_exists = False
    if mpi.is_master_node():
        try:
            os.kill(vasp_pid, 0)
        except OSError as e:
            pid_exists = e.errno == errno.EPERM
        else:
            pid_exists = True

    pid_exists = mpi.bcast(pid_exists)
    return pid_exists

def get_dft_energy():
    """
    Reads energy from the last line of OSZICAR.
    """
    with open('OSZICAR', 'r') as f:
        nextline = f.readline()
        while nextline.strip():
            line = nextline
            nextline = f.readline()
#            print "OSZICAR: ", line[:-1]

    try:
        dft_energy = float(line.split()[2])
    except ValueError:
        print("Cannot read energy from OSZICAR, setting it to zero")
        dft_energy = 0.0

    return dft_energy

class bcolors:
    MAGENTA = '\033[95m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    ENDC = '\033[0m'
#
# Main self-consistent cycle
#
def run_all(vasp_pid, dmft_cycle, cfg_file, n_iter, n_iter_dft, vasp_version):
    """
    """
    mpi.report("  Waiting for VASP lock to appear...")
    while not is_vasp_lock_present():
        time.sleep(1)

    vasp_running = True

    iter = 0
    while vasp_running:
        if debug: print(bcolors.RED + "rank %s"%(mpi.rank) + bcolors.ENDC)
        mpi.report("  Waiting for VASP lock to disappear...")
        mpi.barrier()
        while is_vasp_lock_present():
            time.sleep(1)
#            if debug: print bcolors.YELLOW + " waiting: rank %s"%(mpi.rank) + bcolors.ENDC
            if not is_vasp_running(vasp_pid):
                mpi.report("  VASP stopped")
                vasp_running = False
                break
            
# Tell VASP to stop if the maximum number of iterations is reached
        
        
        if debug: print(bcolors.MAGENTA + "rank %s"%(mpi.rank) + bcolors.ENDC)
        err = 0
        exc = None
        if debug: print(bcolors.BLUE + "plovasp: rank %s"%(mpi.rank) + bcolors.ENDC)
        if mpi.is_master_node():
            converter.generate_and_output_as_text(cfg_file, vasp_dir='./')
            # Read energy from OSZICAR
            dft_energy = get_dft_energy()
        mpi.barrier()

        if debug: print(bcolors.GREEN + "rank %s"%(mpi.rank) + bcolors.ENDC)
        corr_energy, dft_dc = dmft_cycle()
        mpi.barrier()

        if mpi.is_master_node():
            total_energy = dft_energy + corr_energy - dft_dc
            print()
            print("="*80)
            print("  Total energy: ", total_energy)
            print("  DFT energy: ", dft_energy)
            print("  Corr. energy: ", corr_energy)
            print("  DFT DC: ", dft_dc)
            print("="*80)
            print()
        
        # check if we should do additional VASP calculations
        # in the standard VASP version, VASP writes out GAMMA itself
        # so that if we want to keep GAMMA fixed we have to copy it to
        # GAMMA_recent and copy it back after VASP has completed an iteration
        # if we are using a hacked Version of VASP where the write out is
        # disabled we can skip this step.
        # the hack consists of removing the call of LPRJ_LDApU in VASP src file
        # electron.F around line 644
        iter_dft = 0
        
        if vasp_version == 'standard':
            copyfile(src='GAMMA',dst='GAMMA_recent')
        while iter_dft < n_iter_dft:
            if mpi.is_master_node():
                open('./vasp.lock', 'a').close()
            while is_vasp_lock_present():
                time.sleep(1)
                if not is_vasp_running(vasp_pid):
                    mpi.report("  VASP stopped")
                    vasp_running = False
                    break
            iter_dft += 1
            if vasp_version == 'standard':
                copyfile(src='GAMMA_recent',dst='GAMMA')
        iter += 1
        if iter == n_iter:
            print("\n  Maximum number of iterations reached.")
            print("  Aborting VASP iterations...\n")
            f_stop = open('STOPCAR', 'wt')
            f_stop.write("LABORT = .TRUE.\n")
            f_stop.close()
    if mpi.is_master_node():
        total_energy = dft_energy + corr_energy - dft_dc
        with open('TOTENERGY', 'w') as f:
            f.write("  Total energy: %s\n"%(total_energy))
            f.write("  DFT energy: %s\n"%(dft_energy))
            f.write("  Corr. energy: %s\n"%(corr_energy))
            f.write("  DFT DC: %s\n"%(dft_dc))
            f.write("  Energy correction: %s\n"%(corr_energy - dft_dc))

    mpi.report("***Done")

def main():
    
    import importlib
    
    try:
        vasp_pid = int(sys.argv[1])
    except (ValueError, KeyError):
        if mpi.is_master_node():
            print("ERROR: VASP process pid must be provided as the first argument")
        raise

    try:
        n_iter = int(sys.argv[2])
    except (ValueError, KeyError):
        if mpi.is_master_node():
            print("ERROR: Number of iterations must be provided as the second argument")
        raise
        
    try:
        n_iter_dft = int(sys.argv[3])
    except (ValueError, KeyError):
        if mpi.is_master_node():
            print("ERROR: Number of VASP iterations with fixed charge density must be provided as the third argument")
        raise

    try:
        dmft_script = re.sub("\.py$", "", sys.argv[4])
    except:
        if mpi.is_master_node():
            print("ERROR: User-defined DMFT script must be provided as the fourth argument")
        raise

# Optional parameter: config-file name
    try:
        cfg_file = sys.argv[5]
    except KeyError:
        cfg_file = 'plo.cfg'
        
    try:
        vasp_version = sys.argv[6]
    except KeyError:
        vasp_version = 'standard'
    
    if vasp_version != 'standard' and vasp_version != 'no_gamma_write':
        raise Exception('vasp_version has to be standard or no_gamma_write')

#    if len(sys.argv) > 1:
#        vasp_path = sys.argv[1]
#    else:
#        try:
#            vasp_path = os.environ['VASP_DIR']
#        except KeyError:
#            if mpi.is_master_node():
#                print "Path to VASP must be specified either as an argument of in VASP_DIR"
#            raise
    signal.signal(signal.SIGINT, sigint_handler)

    dmft_mod = importlib.import_module(dmft_script)

    run_all(vasp_pid, dmft_mod.dmft_cycle, cfg_file, n_iter, n_iter_dft, vasp_version)

if __name__ == '__main__':
    main()
