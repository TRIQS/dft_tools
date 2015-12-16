r"""
    vasp.main
    =========

    Main script of PLOVasp.

    Runs routines in proper order to generate and store PLOs.

    Usage: python main.py <conf-file> [<path-to-vasp-calcultaion>]
"""
import sys
import vaspio
from inpconf import ConfigParameters
from elstruct import ElectronicStructure
from plotools import generate_plo, output_as_text

def generate_and_output_as_text(conf_filename, vasp_dir):
    """
    Parse config file, process VASP data, and store as text.
    """
# Prepare input-file parameters
    pars = ConfigParameters(conf_filename, verbosity=0)
    pars.parse_input()

# Read VASP data
    if 'efermi' in pars.general:
        efermi_required = False
    else:
        efermi_required = True
    vasp_data = vaspio.VaspData(vasp_dir, efermi_required=efermi_required)
    el_struct = ElectronicStructure(vasp_data)
    el_struct.debug_density_matrix()
    if 'efermi' in pars.general:
        el_struct.efermi = pars.general['efermi']

# Generate and store PLOs
    pshells, pgroups = generate_plo(pars, el_struct)
    output_as_text(pars, el_struct, pshells, pgroups)

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

    generate_and_output_as_text(conf_filename, vasp_dir)

