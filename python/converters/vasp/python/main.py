
import sys
import vaspio
from inpconf import ConfigParameters
from elstruct import ElectronicStructure
from plotools import generate_ortho_plos, plo_output

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


    pars = ConfigParameters(filename, verbosity=0)
    pars.parse_input()
    vasp_data = vaspio.VaspData(vasp_dir)
    el_struct = ElectronicStructure(vasp_data)
    pshells, pgroups = generate_plo(pars, el_struct)
    for gr in pgroups:
        gr.orthogonalize()
    plo_output(pshells, pgroups, el_struct)
