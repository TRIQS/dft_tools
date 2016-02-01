
import os
import errno
import signal
import sys
import time
import pytriqs.utility.mpi as mpi
import converters.plovasp.main as plovasp
from test_ham_hf import dmft_cycle

debug = True
#
# Helper functions
#
def sigint_handler(signal, frame):
    raise SystemExit(1)

def is_vasp_lock_present():
    return os.path.isfile('./vasp.lock')

def is_vasp_running(vasp_pid):
    """
    Tests if VASP initial process is still alive.
    """
    pid_exists = False
    if mpi.is_master_node():
        try:
            os.kill(vasp_pid, 0)
        except OSError, e:
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
        print "Cannot read energy from OSZICAR, setting it to zero"
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
def run_all(vasp_pid):
    """
    """
    mpi.report("  Waiting for VASP lock to appear...")
    while not is_vasp_lock_present():
        time.sleep(1)

    vasp_running = True

    while vasp_running:
        if debug: print bcolors.RED + "rank %s"%(mpi.rank) + bcolors.ENDC
        mpi.report("  Waiting for VASP lock to disappear...")
        mpi.barrier()
        while is_vasp_lock_present():
            time.sleep(1)
#            if debug: print bcolors.YELLOW + " waiting: rank %s"%(mpi.rank) + bcolors.ENDC
            if not is_vasp_running(vasp_pid):
                mpi.report("  VASP stopped")
                vasp_running = False
                break

        if debug: print bcolors.MAGENTA + "rank %s"%(mpi.rank) + bcolors.ENDC
        err = 0
        exc = None
        try:
            if debug: print bcolors.BLUE + "plovasp: rank %s"%(mpi.rank) + bcolors.ENDC
            if mpi.is_master_node():
                plovasp.generate_and_output_as_text('plo.cfg', vasp_dir='./')
# Read energy from OSZICAR
                dft_energy = get_dft_energy()
        except Exception, exc:
            err = 1

        err = mpi.bcast(err)
        if err:
            if mpi.is_master_node():
                raise exc
            else:
                raise SystemExit(1)

        mpi.barrier()

        try:
            if debug: print bcolors.GREEN + "rank %s"%(mpi.rank) + bcolors.ENDC
            corr_energy, dft_dc = dmft_cycle()
        except:
            if mpi.is_master_node():
                print "  master forwarding the exception..."
                raise
            else:
                print "  rank %i exiting..."%(mpi.rank)
                raise SystemExit(1)
        mpi.barrier()

        if mpi.is_master_node():
            total_energy = dft_energy + corr_energy - dft_dc
            print
            print "="*80
            print "  Total energy: ", total_energy
            print "  DFT energy: ", dft_energy
            print "  Corr. energy: ", corr_energy
            print "  DFT DC: ", dft_dc
            print "="*80
            print

        if mpi.is_master_node() and vasp_running:
            open('./vasp.lock', 'a').close()
            
    if mpi.is_master_node():
        total_energy = dft_energy + corr_energy - dft_dc
        with open('TOTENERGY', 'w') as f:
            f.write("  Total energy: %s\n"%(total_energy))
            f.write("  DFT energy: %s\n"%(dft_energy))
            f.write("  Corr. energy: %s\n"%(corr_energy))
            f.write("  DFT DC: %s\n"%(dft_dc))
            f.write("  Energy correction: %s\n"%(corr_energy - dft_dc))

    mpi.report("***Done")


if __name__ == '__main__':
    try:
        vasp_pid = int(sys.argv[1])
    except (ValueError, KeyError):
        if mpi.is_master_node():
            print "VASP process pid must be provided as the first argument"
        raise
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

    run_all(vasp_pid)
