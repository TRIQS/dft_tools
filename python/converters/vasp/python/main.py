
import sys
import vaspio
from inpconf import ConfigParameters

if __name__ == '__main__':
    narg = len(sys.argv)
    if narg < 2:
        raise SystemExit("  Usage: python main.py <conf-file> [<path-to-vasp-calcultaion>]")
    else:
        filename = sys.argv[1]
        if narg > 2:
            vasp_dir = sys.argv[2]
            if vasp_dir[-1] != '/':
                vasp_dir += '/'
        else:
            vasp_dir = './'


    pars = ConfigParameters(filename)
    vasp_data = vaspio.VaspData(vasp_dir)
