
import numpy as np

class ElectronicStructure:
    """
    Class containing electronic structure data.

    **Parameters:**

    - *natom* (int) : total number of atoms
    - *nktot* (int) : total number of `k`-points
    - *nband* (int) : total number of bands
    - *nspin* (int) : spin-polarization
    - *nc_flag* (True/False) : non-collinearity flag
    - *efermi* (float) : Fermi level read from DOSCAR
    - *proj_raw* (array[complex]) : raw projectors from PLOCAR
    - *eigvals* (array[float]) : KS eigenvalues
    - *kmesh* (dict) : parameters of the `k`-mesh
    - *structure* (dict) : parameters of the crystal structure
    - *symmetry* (dict) : paramters of symmetry

    When the object is created a simple consistency check
    of the data coming from different VASP files is performed. 
    """

    def __init__(self, vasp_data):
        self.natom = vasp_data.poscar.nq
        self.nktot = vasp_data.kpoints.nktot

        self.kmesh = {'nktot': self.nktot}
        self.kmesh['kpoints'] = vasp_data.kpoints.kpts
        self.kmesh['kwights'] = vasp_data.eigenval.kwghts
        try:
            self.kmesh['ntet'] = vasp_data.kpoints.ntet
            self.kmesh['itet'] = vasp_data.kpoints.itet
            self.kmesh['volt'] = vasp_data.kpoints.volt
        except NameError:
            pass

# Note that one should not subtract this Fermi level from eigenvalues
# here because the true Fermi level might be provided by conf-file
# (for instance, for spaghetti calculations)
        self.efermi = vasp_data.doscar.efermi

# Note that the number of spin-components of projectors might be different from those
# of bands in case of non-collinear calculations
        self.nspin = vasp_data.plocar.params['ns']
        self.nc_flag = vasp_data.plocar.params['nc_flag'] == 1

        self.nband = vasp_data.eigenval.nband

        self.eigvals = vasp_data.eigenval.eigs

# For later use it is more convenient to use a different order of indices
# [see ProjectorGroup.orthogonalization()]
        self.proj_raw = vasp_data.plocar.plo.transpose((0, 1, 2, 4, 3))

# Check that the number of atoms is the same in PLOCAR and POSCAR
        natom_plo = vasp_data.plocar.params['nion']
        assert natom_plo == self.natom, "PLOCAR is inconsistent with POSCAR (number of atoms)"

# Check that the number of k-points is the same in all files
        nk_plo = vasp_data.plocar.params['nk']
        assert nk_plo == self.nktot, "PLOCAR is inconsistent with IBZKPT (number of k-points)"
        nk_eig = vasp_data.eigenval.nktot
        assert nk_eig == self.nktot, "PLOCAR is inconsistent with EIGENVAL (number of k-points)"

# Check that the number of band is the same in PLOCAR and EIGENVAL
        nb_plo = vasp_data.plocar.params['nb']
        assert nb_plo == self.nband, "PLOCAR is inconsistent with EIGENVAL (number of bands)"



