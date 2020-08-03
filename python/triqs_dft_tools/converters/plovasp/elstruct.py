
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
r"""
    plovasp.elstruct
    ================

    Internal representation of VASP electronic structure data.
"""
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
    - *ferw* (array[float]) : Fermi weights from VASP
    - *kmesh* (dict) : parameters of the `k`-mesh
    - *structure* (dict) : parameters of the crystal structure
    - *symmetry* (dict) : paramters of symmetry

    When the object is created a simple consistency check
    of the data coming from different VASP files is performed.
    """

    def __init__(self, vasp_data):
        self.natom = vasp_data.poscar.nq
        self.type_of_ion = vasp_data.poscar.type_of_ion
        self.nktot = vasp_data.kpoints.nktot

        self.kmesh = {'nktot': self.nktot}
        self.kmesh['kpoints'] = vasp_data.kpoints.kpts
        self.kmesh['kweights'] = vasp_data.kpoints.kwghts
        try:
            self.kmesh['ntet'] = vasp_data.kpoints.ntet
            self.kmesh['itet'] = vasp_data.kpoints.itet
            self.kmesh['volt'] = vasp_data.kpoints.volt
        except AttributeError:
            pass

# Note that one should not subtract this Fermi level from eigenvalues
# here because the true Fermi level might be provided by conf-file
# (for instance, for spaghetti calculations)
        try:
            self.efermi = vasp_data.doscar.efermi
        except AttributeError:
            pass

# Note that the number of spin-components of projectors might be different from those
# of bands in case of non-collinear calculations
        self.nspin = vasp_data.plocar.nspin
        self.nc_flag = vasp_data.plocar.ncdij == 4

        self.nband = vasp_data.plocar.nband

# Check that the number of k-points is the same in all files
        _, ns_plo, nk_plo, nb_plo = vasp_data.plocar.plo.shape
        assert nk_plo == self.nktot, "PLOCAR is inconsistent with IBZKPT (number of k-points)"

# FIXME: Reading from EIGENVAL is obsolete and should be
#        removed completely.
#        if not vasp_data.eigenval.eigs is None:
        if False:
            print("eigvals from EIGENVAL")
            self.eigvals = vasp_data.eigenval.eigs
            self.ferw = vasp_data.eigenval.ferw.transpose((2, 0, 1))

            nk_eig = vasp_data.eigenval.nktot
            assert nk_eig == self.nktot, "PLOCAR is inconsistent with EIGENVAL (number of k-points)"

# Check that the number of band is the same in PROJCAR and EIGENVAL
            assert nb_plo == self.nband, "PLOCAR is inconsistent with EIGENVAL (number of bands)"
        else:
            print("eigvals from LOCPROJ")
            self.eigvals = vasp_data.plocar.eigs
            self.ferw = vasp_data.plocar.ferw.transpose((2, 0, 1))
            self.efermi = vasp_data.doscar.efermi

# For later use it is more convenient to use a different order of indices
# [see ProjectorGroup.orthogonalization()]
        self.proj_raw = vasp_data.plocar.plo
        self.proj_params = vasp_data.plocar.proj_params

# Not needed any more since PROJCAR contains projectors only for a subset of sites
# Check that the number of atoms is the same in PLOCAR and POSCAR
#        natom_plo = vasp_data.plocar.params['nion']
#        assert natom_plo == self.natom, "PLOCAR is inconsistent with POSCAR (number of atoms)"
        self.structure = {'a_brav': vasp_data.poscar.a_brav}
        self.structure['nqtot'] = vasp_data.poscar.nq
        self.structure['kpt_basis'] = vasp_data.poscar.kpt_basis
        self.structure['ntypes'] = vasp_data.poscar.ntypes
        self.structure['nq_types'] = vasp_data.poscar.nions
# Concatenate coordinates grouped by type into one array
        self.structure['qcoords'] = np.vstack(vasp_data.poscar.q_types)
        self.structure['type_of_ion'] = vasp_data.poscar.type_of_ion

        self.kmesh['kpoints_cart'] = 0.0 * self.kmesh['kpoints']

        for ik in range(self.nktot):
            for ii in range(3):
                self.kmesh['kpoints_cart'][ik] += self.kmesh['kpoints'][ik,ii]*self.structure['kpt_basis'][:,ii]

# FIXME: This can be removed if ion coordinates are stored in a continuous array
## Construct a map to access coordinates by index
#        self.structure['ion_index'] = []
#        for isort, nq in enumerate(self.structure['nq_types']):
#            for iq in range(nq):
#                self.structure['ion_index'].append((isort, iq))


    def debug_density_matrix(self):
        """
        Calculate and output the density and overlap matrix out of projectors defined in el_struct.
        """
        plo = self.proj_raw
        nproj, ns, nk, nb = plo.shape
        ions = sorted(list(set([param['isite'] for param in self.proj_params])))
        nions = len(ions)
        norb = nproj // nions

# Spin factor
        sp_fac = 2.0 if ns == 1 and not self.nc_flag else 1.0

        den_mat = np.zeros((ns, nproj, nproj), dtype=np.float64)
        overlap = np.zeros((ns, nproj, nproj), dtype=np.float64)
#        ov_min = np.ones((ns, nproj, nproj), dtype=np.float64) * 100.0
#        ov_max = np.zeros((ns, nproj, nproj), dtype=np.float64)
        for ispin in range(ns):
            for ik in range(nk):
                kweight = self.kmesh['kweights'][ik]
                occ = self.ferw[ispin, ik, :]
                den_mat[ispin, :, :] += np.dot(plo[:, ispin, ik, :] * occ, plo[:, ispin, ik, :].T.conj()).real * kweight * sp_fac
                ov = np.dot(plo[:, ispin, ik, :], plo[:, ispin, ik, :].T.conj()).real
                overlap[ispin, :, :] += ov * kweight
#                ov_max = np.maximum(ov, ov_max)
#                ov_min = np.minimum(ov, ov_min)

# Output only the site-diagonal parts of the matrices
        print()
        print("  Unorthonormalized density matrices and overlaps:")
        for ispin in range(ns):
            print("  Spin:", ispin + 1)
            for io, ion in enumerate(ions):
                print("  Site:", ion)
                iorb_inds = [(ip, param['m']) for ip, param in enumerate(self.proj_params) if param['isite'] == ion]
                norb = len(iorb_inds)
                dm = np.zeros((norb, norb))
                ov = np.zeros((norb, norb))
                for ind, iorb in iorb_inds:
                    for ind2, iorb2 in iorb_inds:
                        dm[iorb, iorb2] = den_mat[ispin, ind, ind2]
                        ov[iorb, iorb2] = overlap[ispin, ind, ind2]

                print("  Density matrix" + (12*norb - 12 + 2)*" " + "Overlap")
                for drow, dov in zip(dm, ov):
                    out = ''.join(map("{0:12.7f}".format, drow))
                    out += "    "
                    out += ''.join(map("{0:12.7f}".format, dov))
                    print(out)
