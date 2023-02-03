
##########################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2018 by G. J. Kraberger
# Copyright (C) 2011 by M. Aichhorn, L. Pourovskii, V. Vildosola
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
##########################################################################
"""
General SumK class and helper functions for combining ab-initio code and triqs
"""

from types import *
import numpy as np
import triqs.utility.dichotomy as dichotomy
from triqs.gf import *
import triqs.utility.mpi as mpi
from triqs.utility.comparison_tests import assert_arrays_are_close
from h5 import *
from .symmetry import *
from .block_structure import BlockStructure
from itertools import product
from warnings import warn
from scipy import compress
from scipy.optimize import minimize, newton, brenth

def compute_DC_from_density(N_tot, U, J, N_spin = None,  n_orbitals=5,  method='cFLL'):
    """
    Computes the double counting correction in terms using various methods


    Parameters
    ----------
    N_tot : float 
        Total density of the impurity
    N_spin : float , default = None
        Spin density, defaults to N_tot*0.5 if not specified
    U : float 
        U value
    J : float 
        J value
    n_orbitals : int, default = 5
        Total number of orbitals
    method : string, default = 'cFLL' 
        possibilities:
        -    cFLL: DC potential from Ryee for spin unpolarized DFT: (DOI: 10.1038/s41598-018-27731-4)
        -    sFLL: same as above for spin polarized DFT
        -    cAMF: around mean field
        -    sAMF: spin polarized around mean field
        -    cHeld: unpolarized Held's formula as reported in (DOI: 10.1103/PhysRevResearch.2.03308)
        -    sHeld: NOT IMPLEMENTED
    
    Returns
    -------
    List of floats:
        -   DC_val: double counting potential
        -   E_val: double counting energy
    """
    if N_spin is not None:
        N_spin2 = N_tot-N_spin
        Mag = N_spin - N_spin2
    L_orbit = (n_orbitals-1)/2 

    match method:
        case 'cFLL':
            E_val = 0.5 * U * N_tot * (N_tot-1) - 0.5 * J * N_tot * (N_tot*0.5-1)
            DC_val = U * (N_tot-0.5) - J *(N_tot*0.5-0.5)

        case 'sFLL':
            assert N_spin is not None, "Spin density not given"
            E_val = 0.5 * U * N_tot * (N_tot-1) - 0.5 * J * N_tot * (N_tot*0.5-1) - 0.25 * J * Mag**2
            DC_val = U * (N_tot-0.5) - J *(N_spin-0.5)
        
        case 'cAMF':
            E_val = +0.5 * U * N_tot **2
            E_val -= 0.25*(U+2*L_orbit*J)/(2*L_orbit+1)*N_tot**2
            DC_val = U * N_tot - 0.5*(U+2*L_orbit*J)/(2*L_orbit+1)*N_tot
            

        case 'sAMF':
            assert N_spin is not None, "Spin density not given"
            E_val = 0.5 * U * N_tot **2
            E_val -= 0.25*(U+2*L_orbit*J)/(2*L_orbit+1)*N_tot**2
            E_val -= 0.25*(U+2*L_orbit*J)/(2*L_orbit+1)*Mag**2
            DC_val = U * N_tot  - (U+2*L_orbit*J)/(2*L_orbit+1)*N_spin
        
        case 'cHeld':
            U_mean = (U + (n_orbitals-1)*(U-2*J)+(n_orbitals-1)*(U-3*J))/(2*n_orbitals-1)
            E_val = 0.5 * U_mean * N_tot * (N_tot - 1)
            DC_val = U_mean * (N_tot-0.5)
        
        case 'sHeld':
            raise ValueError(f"Method sHeld not yet implemented")
        
        case _:
            raise ValueError(f"DC type {method} not supported")

    mpi.report(f"DC potential computed using the {method} method, V_DC = {DC_val:.6f} eV")
    mpi.report(f"E_DC using the {method} method, E_DC = {E_val:.6f} eV")
    if 'Held' in method:
        mpi.report(f"Held method for {n_orbitals} orbitals, computed U_mean={U_mean:.6f} eV")

    return DC_val, E_val


class SumkDFT(object):
    """This class provides a general SumK method for combining ab-initio code and triqs."""

    def __init__(self, hdf_file, h_field=0.0, mesh=None, beta=40, n_iw=1025, use_dft_blocks=False,
                 dft_data='dft_input', symmcorr_data='dft_symmcorr_input', parproj_data='dft_parproj_input',
                 symmpar_data='dft_symmpar_input', bands_data='dft_bands_input', transp_data='dft_transp_input',
                 misc_data='dft_misc_input',bc_data='dft_bandchar_input',fs_data='dft_fs_input'):
        r"""
        Initialises the class from data previously stored into an hdf5 archive.

        Parameters
        ----------
        hdf_file : string
                   Name of hdf5 containing the data.
        h_field : scalar, optional
                  The value of magnetic field to add to the DFT Hamiltonian.
                  The contribution -h_field*sigma is added to diagonal elements of the Hamiltonian.
                  It cannot be used with the spin-orbit coupling on; namely h_field is set to 0 if self.SO=True.
        mesh: MeshImFreq or MeshReFreq, optional. Frequency mesh of Sigma.
        beta : real, optional
               Inverse temperature. Used to construct imaginary frequency if mesh is not given.
        n_iw : integer, optional
               Number of Matsubara frequencies. Used to construct imaginary frequency if mesh is not given.
        use_dft_blocks : boolean, optional
                         If True, the local Green's function matrix for each spin is divided into smaller blocks
                          with the block structure determined from the DFT density matrix of the corresponding correlated shell.

                         Alternatively and additionally, the block structure can be analysed using :meth:`analyse_block_structure <dft.sumk_dft.SumkDFT.analyse_block_structure>`
                         and manipulated using the SumkDFT.block_structre attribute (see :class:`BlockStructure <dft.block_structure.BlockStructure>`).
        dft_data : string, optional
                   Name of hdf5 subgroup in which DFT data for projector and lattice Green's function construction are stored.
        symmcorr_data : string, optional
                        Name of hdf5 subgroup in which DFT data on symmetries of correlated shells
                        (symmetry operations, permutaion matrices etc.) are stored.
        parproj_data : string, optional
                       Name of hdf5 subgroup in which DFT data on non-normalized projectors for non-correlated
                       states (used in the partial density of states calculations) are stored.
        symmpar_data : string, optional
                       Name of hdf5 subgroup in which DFT data on symmetries of the non-normalized projectors
                       are stored.
        bands_data : string, optional
                     Name of hdf5 subgroup in which DFT data necessary for band-structure/k-resolved spectral
                     function calculations (projectors, DFT Hamiltonian for a chosen path in the Brillouin zone etc.)
                     are stored.
        transp_data : string, optional
                      Name of hdf5 subgroup in which DFT data necessary for transport calculations are stored.
        misc_data : string, optional
                    Name of hdf5 subgroup in which miscellaneous DFT data are stored.
        """

        if not isinstance(hdf_file, str):
            mpi.report("Give a string for the hdf5 filename to read the input!")
        else:
            self.hdf_file = hdf_file
            self.dft_data = dft_data
            self.symmcorr_data = symmcorr_data
            self.parproj_data = parproj_data
            self.symmpar_data = symmpar_data
            self.bands_data = bands_data
            self.transp_data = transp_data
            self.misc_data = misc_data
            self.bc_data = bc_data
            self.fs_data = fs_data
            self.h_field = h_field

            if mesh is None:
                self.mesh = MeshImFreq(beta=beta, S='Fermion', n_max=n_iw)
                self.mesh_values = np.linspace(self.mesh(self.mesh.first_index()),
                                               self.mesh(self.mesh.last_index()),
                                               len(self.mesh))
            elif isinstance(mesh, MeshImFreq):
                self.mesh = mesh
                self.mesh_values = np.linspace(self.mesh(self.mesh.first_index()),
                                               self.mesh(self.mesh.last_index()),
                                               len(self.mesh))
            elif isinstance(mesh, MeshReFreq):
                self.mesh = mesh
                self.mesh_values = np.linspace(self.mesh.omega_min, self.mesh.omega_max, len(self.mesh))
            else:
                raise ValueError('mesh must be a triqs mesh of type MeshImFreq or MeshReFreq')

            self.block_structure = BlockStructure()

            # Read input from HDF:
            req_things_to_read = ['energy_unit', 'n_k', 'k_dep_projection', 'SP', 'SO', 'charge_below', 'density_required',
                              'symm_op', 'n_shells', 'shells', 'n_corr_shells', 'corr_shells', 'use_rotations', 'rot_mat',
                              'rot_mat_time_inv', 'n_reps', 'dim_reps', 'T', 'n_orbitals', 'proj_mat', 'bz_weights', 'hopping',
                              'n_inequiv_shells', 'corr_to_inequiv', 'inequiv_to_corr']
            self.subgroup_present, self.values_not_read = self.read_input_from_hdf(
                subgrp=self.dft_data, things_to_read=req_things_to_read)

            # test if all required properties have been found
            if len(self.values_not_read) > 0 and mpi.is_master_node:
                raise ValueError('ERROR: One or more necessary SumK input properties have not been found in the given h5 archive:', self.values_not_read)

            # optional properties to load
            # soon bz_weights is depraced and replaced by kpt_weights, kpts_basis and kpts will become required to read soon
            optional_things_to_read = ['proj_mat_csc', 'proj_or_hk', 'kpt_basis', 'kpts', 'kpt_weights', 'dft_code']
            subgroup_present, self.optional_values_not_read = self.read_input_from_hdf(subgrp=self.dft_data, things_to_read=optional_things_to_read)

            # warning if dft_code was not read (old h5 structure)
            if 'dft_code' in self.optional_values_not_read:
                self.dft_code = None
                mpi.report('\nWarning: old h5 archive without dft_code input flag detected. Please specify sumk.dft_code manually!\n')

            if self.symm_op:
                self.symmcorr = Symmetry(hdf_file, subgroup=self.symmcorr_data)

            if self.SO and (abs(self.h_field) > 0.000001):
                self.h_field = 0.0
                mpi.report(
                    "For SO, the external magnetic field is not implemented, setting it to 0!")

            self.spin_block_names = [['up', 'down'], ['ud']]
            self.n_spin_blocks = [2, 1]
            # Convert spin_block_names to indices -- if spin polarized,
            # differentiate up and down blocks
            self.spin_names_to_ind = [{}, {}]
            for iso in range(2):  # SO = 0 or 1
                for isp in range(self.n_spin_blocks[iso]):
                    self.spin_names_to_ind[iso][
                        self.spin_block_names[iso][isp]] = isp * self.SP

            # GF structure used for the local things in the k sums
            # Most general form allowing for all hybridisation, i.e. largest
            # blocks possible
            self.gf_struct_sumk = [[(sp, self.corr_shells[icrsh]['dim']) for sp in self.spin_block_names[self.corr_shells[icrsh]['SO']]]
                                   for icrsh in range(self.n_corr_shells)]
            # First set a standard gf_struct solver:
            self.gf_struct_solver = [dict([(sp, self.corr_shells[self.inequiv_to_corr[ish]]['dim'])
                                           for sp in self.spin_block_names[self.corr_shells[self.inequiv_to_corr[ish]]['SO']]])
                                     for ish in range(self.n_inequiv_shells)]
            # Set standard (identity) maps from gf_struct_sumk <->
            # gf_struct_solver
            self.sumk_to_solver = [{} for ish in range(self.n_inequiv_shells)]
            self.solver_to_sumk = [{} for ish in range(self.n_inequiv_shells)]
            self.solver_to_sumk_block = [{}
                                         for ish in range(self.n_inequiv_shells)]
            for ish in range(self.n_inequiv_shells):
                for block, inner_dim in self.gf_struct_sumk[self.inequiv_to_corr[ish]]:
                    self.solver_to_sumk_block[ish][block] = block
                    for inner in range(inner_dim):
                        self.sumk_to_solver[ish][
                            (block, inner)] = (block, inner)
                        self.solver_to_sumk[ish][
                            (block, inner)] = (block, inner)
            # assume no shells are degenerate
            self.deg_shells = [[] for ish in range(self.n_inequiv_shells)]

            self.chemical_potential = 0.0  # initialise mu
            self.init_dc()  # initialise the double counting

            # charge mixing parameters
            self.charge_mixing = False
            # defaults from PRB 90 235103 ("... slow but stable convergence ...")
            self.charge_mixing_alpha = 0.1
            self.charge_mixing_gamma = 1.0
            self.deltaNOld = None

            # Analyse the block structure and determine the smallest gf_struct
            # blocks and maps, if desired
            if use_dft_blocks:
                self.analyse_block_structure()

            self.min_band_energy = None
            self.max_band_energy = None

################
# hdf5 FUNCTIONS
################

    def read_input_from_hdf(self, subgrp, things_to_read):
        r"""
        Reads data from the HDF file. Prints a warning if a requested dataset is not found.

        Parameters
        ----------
        subgrp : string
                 Name of hdf5 file subgroup from which the data are to be read.
        things_to_read : list of strings
                         List of datasets to be read from the hdf5 file.

        Returns
        -------
        subgroup_present : boolean
                           Is the subgrp is present in hdf5 file?
        values_not_read : list of strings
                           List of things that could not be read

        """

        values_not_read = []
        # initialise variables on all nodes to ensure mpi broadcast works at
        # the end
        for it in things_to_read:
            setattr(self, it, None)
        subgroup_present = 0

        if mpi.is_master_node():
            with HDFArchive(self.hdf_file, 'r') as ar:
                if subgrp in ar:
                    subgroup_present = True
                    # first read the necessary things:
                    for it in things_to_read:
                        if it in ar[subgrp]:
                            setattr(self, it, ar[subgrp][it])
                        else:
                            values_not_read.append(it)
                else:
                    if (len(things_to_read) != 0):
                        mpi.report(
                            "Loading failed: No %s subgroup in hdf5!" % subgrp)
                    subgroup_present = False
                    values_not_read = things_to_read

        # now do the broadcasting:
        for it in things_to_read:
            setattr(self, it, mpi.bcast(getattr(self, it)))
        subgroup_present = mpi.bcast(subgroup_present)
        values_not_read = mpi.bcast(values_not_read)

        return subgroup_present, values_not_read

    def save(self, things_to_save, subgrp='user_data'):
        r"""
        Saves data from a list into the HDF file. Prints a warning if a requested data is not found in SumkDFT object.

        Parameters
        ----------
        things_to_save : list of strings
                         List of datasets to be saved into the hdf5 file.
        subgrp : string, optional
                 Name of hdf5 file subgroup in which the data are to be stored.
        """

        if not (mpi.is_master_node()):
            return  # do nothing on nodes
        with HDFArchive(self.hdf_file, 'a') as ar:
            if not subgrp in ar: ar.create_group(subgrp)
            for it in things_to_save:
                if it in [ "gf_struct_sumk", "gf_struct_solver",
                        "solver_to_sumk", "sumk_to_solver", "solver_to_sumk_block"]:
                    warn("It is not recommended to save '{}' individually. Save 'block_structure' instead.".format(it))
                try:
                    ar[subgrp][it] = getattr(self, it)
                except:
                    mpi.report("%s not found, and so not saved." % it)

    def load(self, things_to_load, subgrp='user_data'):
        r"""
        Loads user data from the HDF file. Raises an exeption if a requested dataset is not found.

        Parameters
        ----------
        things_to_read : list of strings
                         List of datasets to be read from the hdf5 file.
        subgrp : string, optional
                 Name of hdf5 file subgroup from which the data are to be read.

        Returns
        -------
        list_to_return : list
                         A list containing data read from hdf5.
        """

        if not (mpi.is_master_node()):
            return  # do nothing on nodes
        with HDFArchive(self.hdf_file, 'r') as ar:
            if not subgrp in ar:
                mpi.report("Loading %s failed!" % subgrp)
            list_to_return = []
            for it in things_to_load:
                try:
                    list_to_return.append(ar[subgrp][it])
                except:
                    raise ValueError("load: %s not found, and so not loaded." % it)
        return list_to_return

################
# CORE FUNCTIONS
################

    def downfold(self, ik, ish, bname, gf_to_downfold, gf_inp, shells='corr', ir=None):
        r"""
        Downfolds a block of the Green's function for a given shell and k-point using the corresponding projector matrices.

        Parameters
        ----------
        ik : integer
             k-point index for which the downfolding is to be done.
        ish : integer
              Shell index of GF to be downfolded.

              - if shells='corr': ish labels all correlated shells (equivalent or not)
              - if shells='all': ish labels only representative (inequivalent) non-correlated shells

        bname : string
                Block name of the target block of the lattice Green's function.
        gf_to_downfold : Gf
                       Block of the Green's function that is to be downfolded.
        gf_inp : Gf
                 FIXME
        shells : string, optional

                 - if shells='corr': orthonormalized projectors for correlated shells are used for the downfolding.
                 - if shells='all': non-normalized projectors for all included shells are used for the downfolding.
                 - if shells='csc': orthonormalized projectors for all shells are used for the downfolding. Used for H(k).

        ir : integer, optional
             Index of equivalent site in the non-correlated shell 'ish', only used if shells='all'.

        Returns
        -------
        gf_downfolded : Gf
                      Downfolded block of the lattice Green's function.
        """

        gf_downfolded = gf_inp.copy()
        # get spin index for proj. matrices
        isp = self.spin_names_to_ind[self.SO][bname]
        n_orb = self.n_orbitals[ik, isp]
        if shells == 'corr':
            dim = self.corr_shells[ish]['dim']
            projmat = self.proj_mat[ik, isp, ish, 0:dim, 0:n_orb]
        elif shells == 'all':
            if ir is None:
                raise ValueError("downfold: provide ir if treating all shells.")
            dim = self.shells[ish]['dim']
            projmat = self.proj_mat_all[ik, isp, ish, ir, 0:dim, 0:n_orb]
        elif shells == 'csc':
            projmat = self.proj_mat_csc[ik, isp, :, 0:n_orb]

        gf_downfolded.from_L_G_R(
            projmat, gf_to_downfold, projmat.conjugate().transpose())

        return gf_downfolded

    def upfold(self, ik, ish, bname, gf_to_upfold, gf_inp, shells='corr', ir=None):
        r"""
        Upfolds a block of the Green's function for a given shell and k-point using the corresponding projector matrices.

        Parameters
        ----------
        ik : integer
             k-point index for which the upfolding is to be done.
        ish : integer
              Shell index of GF to be upfolded.

              - if shells='corr': ish labels all correlated shells (equivalent or not)
              - if shells='all': ish labels only representative (inequivalent) non-correlated shells

        bname : string
                Block name of the target block of the lattice Green's function.
        gf_to_upfold : Gf
                       Block of the Green's function that is to be upfolded.
        gf_inp : Gf
                 FIXME
        shells : string, optional

                 - if shells='corr': orthonormalized projectors for correlated shells are used for the upfolding.
                 - if shells='all': non-normalized projectors for all included shells are used for the upfolding.
                 - if shells='csc': orthonormalized projectors for all shells are used for the upfolding. Used for H(k).

        ir : integer, optional
             Index of equivalent site in the non-correlated shell 'ish', only used if shells='all'.

        Returns
        -------
        gf_upfolded : Gf
                      Upfolded block of the lattice Green's function.
        """

        gf_upfolded = gf_inp.copy()
        # get spin index for proj. matrices
        isp = self.spin_names_to_ind[self.SO][bname]
        n_orb = self.n_orbitals[ik, isp]
        if shells == 'corr':
            dim = self.corr_shells[ish]['dim']
            projmat = self.proj_mat[ik, isp, ish, 0:dim, 0:n_orb]
        elif shells == 'all':
            if ir is None:
                raise ValueError("upfold: provide ir if treating all shells.")
            dim = self.shells[ish]['dim']
            projmat = self.proj_mat_all[ik, isp, ish, ir, 0:dim, 0:n_orb]
        elif shells == 'csc':
            projmat = self.proj_mat_csc[ik, isp, 0:n_orb, 0:n_orb]

        gf_upfolded.from_L_G_R(
            projmat.conjugate().transpose(), gf_to_upfold, projmat)

        return gf_upfolded

    def rotloc(self, ish, gf_to_rotate, direction, shells='corr'):
        r"""
        Rotates a block of the local Green's function from the local frame to the global frame and vice versa.

        Parameters
        ----------
        ish : integer
              Shell index of GF to be rotated.

              - if shells='corr': ish labels all correlated shells (equivalent or not)
              - if shells='all': ish labels only representative (inequivalent) non-correlated shells

        gf_to_rotate : Gf
                       Block of the Green's function that is to be rotated.
        direction : string
                    The direction of rotation can be either

                    - 'toLocal' : global -> local transformation,
                    - 'toGlobal' : local -> global transformation.

        shells : string, optional

                 - if shells='corr': the rotation matrix for the correlated shell 'ish' is used,
                 - if shells='all': the rotation matrix for the generic (non-correlated) shell 'ish' is used.

        Returns
        -------
        gf_rotated : Gf
                     Rotated block of the local Green's function.
        """

        assert ((direction == 'toLocal') or (direction == 'toGlobal')
                ), "rotloc: Give direction 'toLocal' or 'toGlobal'."
        gf_rotated = gf_to_rotate.copy()
        if shells == 'corr':
            rot_mat_time_inv = self.rot_mat_time_inv
            rot_mat = self.rot_mat
        elif shells == 'all':
            rot_mat_time_inv = self.rot_mat_all_time_inv
            rot_mat = self.rot_mat_all

        if direction == 'toGlobal':

            if (rot_mat_time_inv[ish] == 1) and self.SO:
                gf_rotated << gf_rotated.transpose()
                gf_rotated.from_L_G_R(rot_mat[ish].conjugate(
                ), gf_rotated, rot_mat[ish].transpose())
            else:
                gf_rotated.from_L_G_R(rot_mat[ish], gf_rotated, rot_mat[
                                      ish].conjugate().transpose())

        elif direction == 'toLocal':

            if (rot_mat_time_inv[ish] == 1) and self.SO:
                gf_rotated << gf_rotated.transpose()
                gf_rotated.from_L_G_R(rot_mat[ish].transpose(
                ), gf_rotated, rot_mat[ish].conjugate())
            else:
                gf_rotated.from_L_G_R(rot_mat[ish].conjugate(
                ).transpose(), gf_rotated, rot_mat[ish])

        return gf_rotated

    def lattice_gf(self, ik, mu=None, broadening=None, mesh=None, with_Sigma=True, with_dc=True):
        r"""
        Calculates the lattice Green function for a given k-point from the DFT Hamiltonian and the self energy.

        Parameters
        ----------
        ik : integer
             k-point index.
        mu : real, optional
             Chemical potential for which the Green's function is to be calculated.
             If not provided, self.chemical_potential is used for mu.
        broadening : real, optional
                     Imaginary shift for the axis along which the real-axis GF is calculated.
                     If not provided, broadening will be set to double of the distance between mesh points in 'mesh'.
        mesh : MeshReFreq or MeshImFreq, optional
                    Mesh to be used if with_Sigma=False. If with Sigma=False and mesh is none then self.mesh is used.
        with_Sigma : boolean, optional
                     If True the GF will be calculated with the self-energy stored in self.Sigmaimp_(w/iw), for real/Matsubara GF, respectively.
                     In this case the mesh is taken from the self.Sigma_imp object.
                     If with_Sigma=True but self.Sigmaimp_(w/iw) is not present, with_Sigma is reset to False.
        with_dc : boolean, optional
                  if True and with_Sigma=True, the dc correction is substracted from the self-energy before it is included into GF.

        Returns
        -------
        G_latt : BlockGf
                 Lattice Green's function.

        """
        if mu is None:
            mu = self.chemical_potential
        ntoi = self.spin_names_to_ind[self.SO]
        spn = self.spin_block_names[self.SO]
        if not hasattr(self, "Sigma_imp"):
            with_Sigma = False
        if broadening is None:
            if mesh is None:
                broadening = 0.01
            else:  # broadening = 2 * \Delta omega, where \Delta omega is the spacing of omega points
                broadening = 2.0 * ((mesh[1] - mesh[0]) / (mesh[2] - 1))

        # Check if G_latt is present
        set_up_G_latt = False                       # Assume not
        if not hasattr(self, "G_latt" ):
            # Need to create G_latt_(i)w
            set_up_G_latt = True
        else:                                       # Check that existing GF is consistent
            G_latt = self.G_latt
            GFsize = [gf.target_shape[0] for bname, gf in G_latt]
            unchangedsize = all([self.n_orbitals[ik, ntoi[spn[isp]]] == GFsize[
                                isp] for isp in range(self.n_spin_blocks[self.SO])])
            if (not mesh is None) or (not unchangedsize):
                set_up_G_latt = True

        # Are we including Sigma?
        if with_Sigma:
            Sigma_imp = self.Sigma_imp
            if with_dc:
                sigma_minus_dc = self.add_dc()
            else:
                sigma_minus_dc = [s for s in Sigma_imp]
            if not mesh is None:
                warn('lattice_gf called with Sigma and given mesh. Mesh will be taken from Sigma.')

            if self.mesh != Sigma_imp[0].mesh:
                mesh = Sigma_imp[0].mesh
                if isinstance(mesh, MeshImFreq):
                    mesh_values = np.linspace(mesh(mesh.first_index()), mesh(mesh.last_index()), len(mesh))
                else:
                    mesh_values = np.linspace(mesh.omega_min, mesh.omega_max, len(mesh))
            else:
                mesh = self.mesh
                mesh_values = self.mesh_values

        elif not mesh is None:
            assert isinstance(mesh, MeshReFreq) or isinstance(mesh, MeshImFreq),  "mesh must be a triqs MeshReFreq or MeshImFreq"
            if isinstance(mesh, MeshImFreq):
                mesh_values = np.linspace(mesh(mesh.first_index()), mesh(mesh.last_index()), len(mesh))
            else:
                mesh_values = np.linspace(mesh.omega_min, mesh.omega_max, len(mesh))
        else:
            mesh = self.mesh
            mesh_values = self.mesh_values

        # Set up G_latt
        if set_up_G_latt:
            block_structure = [
                list(range(self.n_orbitals[ik, ntoi[sp]])) for sp in spn]
            gf_struct = [(spn[isp], block_structure[isp])
                         for isp in range(self.n_spin_blocks[self.SO])]
            block_ind_list = [block for block, inner in gf_struct]
            if isinstance(mesh, MeshImFreq):
                glist = lambda: [Gf(indices=inner, mesh=mesh)
                                 for block, inner in gf_struct]
            else:
                glist = lambda: [Gf(indices=inner, mesh=mesh)
                                 for block, inner in gf_struct]
            G_latt = BlockGf(name_list=block_ind_list,
                             block_list=glist(), make_copies=False)
            G_latt.zero()

        idmat = [np.identity(
            self.n_orbitals[ik, ntoi[sp]], complex) for sp in spn]

        # fill Glatt
        for ibl, (block, gf) in enumerate(G_latt):
            ind = ntoi[spn[ibl]]
            n_orb = self.n_orbitals[ik, ind]
            if isinstance(mesh, MeshImFreq):
                gf.data[:, :, :] = (idmat[ibl] * (mesh_values[:, None, None] + mu + self.h_field*(1-2*ibl))
                                    - self.hopping[ik, ind, 0:n_orb, 0:n_orb])
            else:
                gf.data[:, :, :] = (idmat[ibl] *
                                    (mesh_values[:, None, None] + mu + self.h_field*(1-2*ibl) + 1j*broadening)
                                    - self.hopping[ik, ind, 0:n_orb, 0:n_orb])

            if with_Sigma:
                for icrsh in range(self.n_corr_shells):
                    gf -= self.upfold(ik, icrsh, block,
                                      sigma_minus_dc[icrsh][block], gf)

        G_latt.invert()
        self.G_latt = G_latt

        return G_latt

    def set_Sigma(self, Sigma_imp, transform_to_sumk_blocks=True):
        self.put_Sigma(Sigma_imp, transform_to_sumk_blocks)

    def put_Sigma(self, Sigma_imp, transform_to_sumk_blocks=True):
        r"""
        Insert the impurity self-energies into the sumk_dft class.

        Parameters
        ----------
        Sigma_imp : list of BlockGf (Green's function) objects
            List containing impurity self-energy for all (inequivalent) correlated shells.
            Self-energies for equivalent shells are then automatically set by this function.
            The self-energies can be of the real or imaginary-frequency type.
        transform_to_sumk_blocks : bool, optional
            If True (default), the input Sigma_imp will be transformed to the block structure ``gf_struct_sumk``,
            else it has to be given in ``gf_struct_sumk``.
        """

        if transform_to_sumk_blocks:
            Sigma_imp = self.transform_to_sumk_blocks(Sigma_imp)

        assert isinstance(Sigma_imp, list),\
            "put_Sigma: Sigma_imp has to be a list of Sigmas for the correlated shells, even if it is of length 1!"
        assert len(Sigma_imp) == self.n_corr_shells,\
            "put_Sigma: give exactly one Sigma for each corr. shell!"

        if isinstance(self.mesh, MeshImFreq) and all(isinstance(gf.mesh, MeshImFreq) and isinstance(gf, Gf) and gf.mesh == self.mesh for bname, gf in Sigma_imp[0]):
            # Imaginary frequency Sigma:
            self.Sigma_imp = [self.block_structure.create_gf(ish=icrsh, mesh=Sigma_imp[icrsh].mesh, space='sumk')
                                 for icrsh in range(self.n_corr_shells)]
            SK_Sigma_imp = self.Sigma_imp
        elif isinstance(self.mesh, MeshReFreq) and all(isinstance(gf, Gf) and isinstance(gf.mesh, MeshReFreq) and gf.mesh == self.mesh for bname, gf in Sigma_imp[0]):
            # Real frequency Sigma:
            self.Sigma_imp = [self.block_structure.create_gf(ish=icrsh, mesh=Sigma_imp[icrsh].mesh, gf_function=Gf, space='sumk')
                                for icrsh in range(self.n_corr_shells)]
            SK_Sigma_imp = self.Sigma_imp

        else:
            raise ValueError("put_Sigma: Sigma_imp must have the same mesh as SumKDFT.mesh.")

        # rotation from local to global coordinate system:
        for icrsh in range(self.n_corr_shells):
            for bname, gf in SK_Sigma_imp[icrsh]:
                if self.use_rotations:
                    gf << self.rotloc(icrsh,
                                      Sigma_imp[icrsh][bname],
                                      direction='toGlobal')
                else:
                    gf << Sigma_imp[icrsh][bname]

        #warning if real frequency self energy is within the bounds of the band energies
        if isinstance(self.mesh, MeshReFreq):
            if self.min_band_energy is None or self.max_band_energy is None:
                self.calculate_min_max_band_energies()
            mesh = np.array([i for i in self.mesh.values()])
            if mesh[0] > (self.min_band_energy - self.chemical_potential) or mesh[-1] < (self.max_band_energy - self.chemical_potential):
                warn('The given Sigma is on a mesh which does not cover the band energy range. The Sigma MeshReFreq runs from %f to %f, while the band energy (minus the chemical potential) runs from %f to %f'%(mesh[0], mesh[-1], self.min_band_energy, self.max_band_energy))

    def transform_to_sumk_blocks(self, Sigma_imp, Sigma_out=None):
        r""" transform Sigma from solver to sumk space

        Parameters
        ----------
        Sigma_imp : list of BlockGf (Green's function) objects
            List containing impurity self-energy for all inequivalent correlated shells.
            The self-energies can be of the real or imaginary-frequency type.
        Sigma_out : list of BlockGf
            list of one BlockGf per correlated shell with the block structure
            according to ``gf_struct_sumk``; if None, it will be created
        """

        assert isinstance(Sigma_imp, list),\
            "transform_to_sumk_blocks: Sigma_imp has to be a list of Sigmas for the inequivalent correlated shells, even if it is of length 1!"
        assert len(Sigma_imp) == self.n_inequiv_shells,\
            "transform_to_sumk_blocks: give exactly one Sigma for each inequivalent corr. shell!"

        if Sigma_out is None:
            Sigma_out = [self.block_structure.create_gf(ish=icrsh, mesh=Sigma_imp[self.corr_to_inequiv[icrsh]].mesh, space='sumk')
                         for icrsh in range(self.n_corr_shells)]
        else:
            for icrsh in range(self.n_corr_shells):
                self.block_structure.check_gf(Sigma_out,
                                              ish=icrsh,
                                              space='sumk')

        # transform the CTQMC blocks to the full matrix:
        for icrsh in range(self.n_corr_shells):
            # ish is the index of the inequivalent shell corresponding to icrsh
            ish = self.corr_to_inequiv[icrsh]
            self.block_structure.convert_gf(
                G=Sigma_imp[ish],
                G_struct=None,
                space_from='solver',
                space_to='sumk',
                ish_from=ish,
                ish_to=icrsh,
                G_out=Sigma_out[icrsh])
        return Sigma_out

    def extract_G_loc(self, mu=None, with_Sigma=True, with_dc=True, broadening=None,
                      transform_to_solver_blocks=True, show_warnings=True):
        r"""
        Extracts the local downfolded Green function by the Brillouin-zone integration of the lattice Green's function.

        Parameters
        ----------
        mu : real, optional
            Input chemical potential. If not provided the value of self.chemical_potential is used as mu.
        with_Sigma : boolean, optional
            If True then the local GF is calculated with the self-energy self.Sigma_imp.
        with_dc : boolean, optional
            If True then the double-counting correction is subtracted from the self-energy in calculating the GF.
        broadening : float, optional
            Imaginary shift for the axis along which the real-axis GF is calculated.
            If not provided, broadening will be set to double of the distance between mesh points in 'mesh'.
            Only relevant for real-frequency GF.
        transform_to_solver_blocks : bool, optional
            If True (default), the returned G_loc will be transformed to the block structure ``gf_struct_solver``,
            else it will be in ``gf_struct_sumk``.
        show_warnings : bool, optional
            Displays warning messages during transformation
            (Only effective if transform_to_solver_blocks = True

        Returns
        -------
        G_loc : list of BlockGf (Green's function) objects
            List of the local Green's functions for all (inequivalent) correlated shells,
            rotated into the corresponding local frames.
            If ``transform_to_solver_blocks`` is True, it will be one per inequivalent correlated shell, else one per
            correlated shell.
        """

        if mu is None:
            mu = self.chemical_potential

        if with_Sigma:
            mesh = self.Sigma_imp[0].mesh
            if mesh != self.mesh:
                warn('self.mesh and self.Sigma_imp[0].mesh are differen! Using mesh from Sigma')

        else:
            mesh = self.mesh

        # create G_loc to be returned in sumk space for all correlated shells. Trafo to solver block structure done later
        G_loc = [self.block_structure.create_gf(ish=ish, mesh=mesh, space='sumk') for ish in range(self.n_corr_shells)]

        ikarray = np.array(list(range(self.n_k)))
        for ik in mpi.slice_array(ikarray):
            if isinstance(self.mesh, MeshImFreq):
                G_latt = self.lattice_gf(
                    ik=ik, mu=mu, with_Sigma=with_Sigma, with_dc=with_dc)
            else:
                G_latt = self.lattice_gf(
                    ik=ik, mu=mu, with_Sigma=with_Sigma, with_dc=with_dc, broadening=broadening)
            G_latt *= self.bz_weights[ik]

            for icrsh in range(self.n_corr_shells):
                # init temporary storage
                for bname, gf in G_loc[icrsh]:
                    gf += self.downfold(ik, icrsh, bname, G_latt[bname], gf)

        # Collect data from mpi
        for icrsh in range(self.n_corr_shells):
            G_loc[icrsh] << mpi.all_reduce(
                mpi.world, G_loc[icrsh], lambda x, y: x + y)
        mpi.barrier()

        # G_loc[:] is now the sum over k projected to the local orbitals.
        # here comes the symmetrisation, if needed:
        if self.symm_op != 0:
            G_loc = self.symmcorr.symmetrize(G_loc)

        # G_loc is rotated to the local coordinate system:
        if self.use_rotations:
            for icrsh in range(self.n_corr_shells):
                for bname, gf in G_loc[icrsh]:
                    G_loc[icrsh][bname] << self.rotloc(
                        icrsh, gf, direction='toLocal')

        if transform_to_solver_blocks:
            return self.transform_to_solver_blocks(G_loc, show_warnings=show_warnings)

        return G_loc

    def transform_to_solver_blocks(self, G_loc, G_out=None, show_warnings = True):
        """ transform G_loc from sumk to solver space

        Parameters
        ----------
        G_loc : list of BlockGf
            a list of one BlockGf per correlated shell with a structure
            according to ``gf_struct_sumk``, e.g. as returned by
            :py:meth:`.extract_G_loc` with ``transform_to_solver_blocks=False``.
        G_out : list of BlockGf
            a list of one BlockGf per *inequivalent* correlated shell
            with a structure according to ``gf_struct_solver``.
            The output Green's function (if not given, a new one is
            created)

        Returns
        -------
        G_out
        """

        assert isinstance(G_loc, list), "G_loc must be a list (with elements for each correlated shell)"

        if G_out is None:
            G_out = [self.block_structure.create_gf(ish=ish, mesh=G_loc[self.inequiv_to_corr[ish]].mesh)
                     for ish in range(self.n_inequiv_shells)]
        else:
            for ish in range(self.n_inequiv_shells):
                self.block_structure.check_gf(G_out, ish=ish)

        # transform to CTQMC blocks:
        for ish in range(self.n_inequiv_shells):
            self.block_structure.convert_gf(
                G=G_loc[self.inequiv_to_corr[ish]],
                G_struct=None,
                ish_from=self.inequiv_to_corr[ish],
                ish_to=ish,
                space_from='sumk',
                G_out=G_out[ish],
                show_warnings = show_warnings)

        # return only the inequivalent shells:
        return G_out

    def analyse_block_structure(self, threshold=0.00001, include_shells=None, dm=None, hloc=None):
        r"""
        Determines the block structure of local Green's functions by analysing the structure of
        the corresponding density matrices and the local Hamiltonian. The resulting block structures
        for correlated shells are stored in the :class:`SumkDFT.block_structure <dft.block_structure.BlockStructure>` attribute.

        Parameters
        ----------
        threshold : real, optional
                    If the difference between density matrix / hloc elements is below threshold,
                    they are considered to be equal.
        include_shells : list of integers, optional
                         List of correlated shells to be analysed.
                         If include_shells is not provided all correlated shells will be analysed.
        dm : list of dict, optional
             List of density matrices from which block stuctures are to be analysed.
             Each density matrix is a dict {block names: 2d numpy arrays} for each correlated shell.
             If not provided, dm will be calculated from the DFT Hamiltonian by a simple-point BZ integration.
        hloc : list of dict, optional
               List of local Hamiltonian matrices from which block stuctures are to be analysed
               Each Hamiltonian is a dict {block names: 2d numpy arrays} for each inequivalent shell.
               If not provided, it will be calculated using eff_atomic_levels.
        """

        self.gf_struct_solver = [{} for ish in range(self.n_inequiv_shells)]
        self.sumk_to_solver = [{} for ish in range(self.n_inequiv_shells)]
        self.solver_to_sumk = [{} for ish in range(self.n_inequiv_shells)]
        self.solver_to_sumk_block = [{}
                                     for ish in range(self.n_inequiv_shells)]

        if dm is None:
            dm = self.density_matrix(method='using_point_integration')
        dens_mat = [dm[self.inequiv_to_corr[ish]]
                    for ish in range(self.n_inequiv_shells)]
        if hloc is None:
            hloc = self.eff_atomic_levels()

        if include_shells is None:
            include_shells = list(range(self.n_inequiv_shells))
        for ish in include_shells:

            for sp in self.spin_block_names[self.corr_shells[self.inequiv_to_corr[ish]]['SO']]:
                n_orb = self.corr_shells[self.inequiv_to_corr[ish]]['dim']
                # gives an index list of entries larger that threshold
                dmbool = (abs(dens_mat[ish][sp]) > threshold)
                hlocbool = (abs(hloc[ish][sp]) > threshold)

                # Determine off-diagonal entries in upper triangular part of
                # density matrix
                offdiag = set([])
                for i in range(n_orb):
                    for j in range(i + 1, n_orb):
                        if dmbool[i, j] or hlocbool[i, j]:
                            offdiag.add((i, j))

                # Determine the number of non-hybridising blocks in the gf
                blocs = [[i] for i in range(n_orb)]
                while len(offdiag) != 0:
                    pair = offdiag.pop()
                    for b1, b2 in product(blocs, blocs):
                        if (pair[0] in b1) and (pair[1] in b2):
                            if blocs.index(b1) != blocs.index(b2):     # In separate blocks?
                                # Merge two blocks
                                b1.extend(blocs.pop(blocs.index(b2)))
                                break                                  # Move on to next pair in offdiag

                # Set the gf_struct for the solver accordingly
                num_blocs = len(blocs)
                for i in range(num_blocs):
                    blocs[i].sort()
                    self.gf_struct_solver[ish].update(
                        [('%s_%s' % (sp, i), len(blocs[i]))])

                # Construct sumk_to_solver taking (sumk_block, sumk_index) --> (solver_block, solver_inner)
                # and solver_to_sumk taking (solver_block, solver_inner) -->
                # (sumk_block, sumk_index)
                for i in range(num_blocs):
                    for j in range(len(blocs[i])):
                        block_sumk = sp
                        inner_sumk = blocs[i][j]
                        block_solv = '%s_%s' % (sp, i)
                        inner_solv = j
                        self.sumk_to_solver[ish][(block_sumk, inner_sumk)] = (
                            block_solv, inner_solv)
                        self.solver_to_sumk[ish][(block_solv, inner_solv)] = (
                            block_sumk, inner_sumk)
                        self.solver_to_sumk_block[ish][block_solv] = block_sumk

            # Now calculate degeneracies of orbitals
            dm = {}
            for block, block_dim in self.gf_struct_solver[ish].items():
                # get dm for the blocks:
                dm[block] = np.zeros(
                    [block_dim, block_dim], complex)
                for ind1 in range(block_dim):
                    for ind2 in range(block_dim):
                        block_sumk, ind1_sumk = self.solver_to_sumk[
                            ish][(block, ind1)]
                        block_sumk, ind2_sumk = self.solver_to_sumk[
                            ish][(block, ind2)]
                        dm[block][ind1, ind2] = dens_mat[ish][
                            block_sumk][ind1_sumk, ind2_sumk]

            for block1 in self.gf_struct_solver[ish].keys():
                for block2 in self.gf_struct_solver[ish].keys():
                    if dm[block1].shape == dm[block2].shape:
                        if ((abs(dm[block1] - dm[block2]) < threshold).all()) and (block1 != block2):
                            ind1 = -1
                            ind2 = -2
                            # check if it was already there:
                            for n, ind in enumerate(self.deg_shells[ish]):
                                if block1 in ind:
                                    ind1 = n
                                if block2 in ind:
                                    ind2 = n
                            if (ind1 < 0) and (ind2 >= 0):
                                self.deg_shells[ish][ind2].append(block1)
                            elif (ind1 >= 0) and (ind2 < 0):
                                self.deg_shells[ish][ind1].append(block2)
                            elif (ind1 < 0) and (ind2 < 0):
                                self.deg_shells[ish].append([block1, block2])

    def _get_hermitian_quantity_from_gf(self, G):
        """ Convert G to a Hermitian quantity

        For G(tau) and G(iw), G(tau) is returned.
        For G(t) and G(w), the spectral function is returned.

        Parameters
        ----------
        G : list of BlockGf of GfImFreq, GfImTime, GfReFreq or GfReTime
            the input Green's function

        Returns
        -------
        gf : list of BlockGf of GfImTime or GfReFreq
            the output G(tau) or A(w)
        """
        # make a GfImTime from the supplied GfImFreq
        if all(isinstance(g_sh.mesh, MeshImFreq) for g_sh in G):
            gf = [BlockGf(name_block_generator = [(name, GfImTime(beta=block.mesh.beta,
                indices=block.indices,n_points=len(block.mesh)+1)) for name, block in g_sh],
                make_copies=False) for g_sh in G]
            for ish in range(len(gf)):
                for name, g in gf[ish]:
                    g.set_from_fourier(G[ish][name])
        # keep a GfImTime from the supplied GfImTime
        elif all(isinstance(g_sh.mesh, MeshImTime) for g_sh in G):
            gf = G
        # make a spectral function from the supplied GfReFreq
        elif all(isinstance(g_sh.mesh, MeshReFreq) for g_sh in G):
            gf = [g_sh.copy() for g_sh in G]
            for ish in range(len(gf)):
                for name, g in gf[ish]:
                    g << 1.0j*(g-g.conjugate().transpose())/2.0/np.pi
        elif all(isinstance(g_sh.mesh, MeshReTime) for g_sh in G):
            def get_delta_from_mesh(mesh):
                w0 = None
                for w in mesh:
                    if w0 is None:
                        w0 = w
                    else:
                        return w-w0
            gf = [BlockGf(name_block_generator = [(name, GfReFreq(
                window=(-np.pi*(len(block.mesh)-1) / (len(block.mesh)*get_delta_from_mesh(block.mesh)),
                np.pi*(len(block.mesh)-1) / (len(block.mesh)*get_delta_from_mesh(block.mesh))),
                n_points=len(block.mesh), indices=block.indices)) for name, block in g_sh], make_copies=False)
                for g_sh in G]

            for ish in range(len(gf)):
                for name, g in gf[ish]:
                    g.set_from_fourier(G[ish][name])
                    g << 1.0j*(g-g.conjugate().transpose())/2.0/np.pi
        else:
            raise Exception("G must be a list of BlockGf of either GfImFreq, GfImTime, GfReFreq or GfReTime")
        return gf



    def analyse_block_structure_from_gf(self, G, threshold=1.e-5, include_shells=None, analyse_deg_shells = True):
        r"""
        Determines the block structure of local Green's functions by analysing
        the structure of the corresponding non-interacting Green's function.
        The resulting block structures for correlated shells are
        stored in the :class:`SumkDFT.block_structure <dft.block_structure.BlockStructure>`
        attribute.

        This is a safer alternative to analyse_block_structure, because
        the full non-interacting Green's function is taken into account
        and not just the density matrix and Hloc.

        Parameters
        ----------
        G : list of BlockGf of GfImFreq, GfImTime, GfReFreq or GfReTime
            the non-interacting Green's function for each inequivalent correlated shell
        threshold : real, optional
                    If the difference between matrix elements is below threshold,
                    they are considered to be equal.
        include_shells : list of integers, optional
                         List of correlated shells to be analysed.
                         If include_shells is not provided all correlated shells will be analysed.
        analyse_deg_shells : bool
                             Whether to call the analyse_deg_shells function
                             after having finished the block structure analysis

        Returns
        -------
        G : list of BlockGf of GfImFreq or GfImTime
            the Green's function transformed into the new block structure
        """

        assert isinstance(G, list), "G must be a list (with elements for each correlated shell)"

        gf = self._get_hermitian_quantity_from_gf(G)

        # initialize the variables
        self.gf_struct_solver = [{} for ish in range(self.n_inequiv_shells)]
        self.sumk_to_solver = [{} for ish in range(self.n_inequiv_shells)]
        self.solver_to_sumk = [{} for ish in range(self.n_inequiv_shells)]
        self.solver_to_sumk_block = [{}
                                     for ish in range(self.n_inequiv_shells)]

        # the maximum value of each matrix element of each block and shell
        max_gf = [{name:np.max(np.abs(g.data),0) for name, g in gf[ish]} for ish in range(self.n_inequiv_shells)]

        if include_shells is None:
            # include all shells
            include_shells = list(range(self.n_inequiv_shells))

        for ish in include_shells:
            for sp in self.spin_block_names[self.corr_shells[self.inequiv_to_corr[ish]]['SO']]:
                n_orb = self.corr_shells[self.inequiv_to_corr[ish]]['dim']

                # gives an index list of entries larger that threshold
                maxgf_bool = (abs(max_gf[ish][sp]) > threshold)

                # Determine off-diagonal entries in upper triangular part of the
                # Green's function
                offdiag = set([])
                for i in range(n_orb):
                    for j in range(i + 1, n_orb):
                        if maxgf_bool[i, j]:
                            offdiag.add((i, j))

                # Determine the number of non-hybridising blocks in the gf
                blocs = [[i] for i in range(n_orb)]
                while len(offdiag) != 0:
                    pair = offdiag.pop()
                    for b1, b2 in product(blocs, blocs):
                        if (pair[0] in b1) and (pair[1] in b2):
                            if blocs.index(b1) != blocs.index(b2):     # In separate blocks?
                                # Merge two blocks
                                b1.extend(blocs.pop(blocs.index(b2)))
                                break                                  # Move on to next pair in offdiag

                # Set the gf_struct for the solver accordingly
                num_blocs = len(blocs)
                for i in range(num_blocs):
                    blocs[i].sort()
                    self.gf_struct_solver[ish].update(
                        [('%s_%s' % (sp, i), len(blocs[i]))])

                # Construct sumk_to_solver taking (sumk_block, sumk_index) --> (solver_block, solver_inner)
                # and solver_to_sumk taking (solver_block, solver_inner) -->
                # (sumk_block, sumk_index)
                for i in range(num_blocs):
                    for j in range(len(blocs[i])):
                        block_sumk = sp
                        inner_sumk = blocs[i][j]
                        block_solv = '%s_%s' % (sp, i)
                        inner_solv = j
                        self.sumk_to_solver[ish][(block_sumk, inner_sumk)] = (
                            block_solv, inner_solv)
                        self.solver_to_sumk[ish][(block_solv, inner_solv)] = (
                            block_sumk, inner_sumk)
                        self.solver_to_sumk_block[ish][block_solv] = block_sumk

        # transform G to the new structure
        full_structure = BlockStructure.full_structure(
            [{sp:self.corr_shells[self.inequiv_to_corr[ish]]['dim']
                for sp in self.spin_block_names[self.corr_shells[self.inequiv_to_corr[ish]]['SO']]}
                for ish in range(self.n_inequiv_shells)],self.corr_to_inequiv)
        G_transformed = [
            self.block_structure.convert_gf(G[ish],
                full_structure, ish, mesh=G[ish].mesh.copy(), show_warnings=threshold,
                gf_function=type(G[ish]._first()), space_from='sumk', space_to='solver')
            for ish in range(self.n_inequiv_shells)]

        if analyse_deg_shells:
            self.analyse_deg_shells(G_transformed, threshold, include_shells)
        return G_transformed

    def analyse_deg_shells(self, G, threshold=1.e-5, include_shells=None):
        r"""
        Determines the degenerate shells of local Green's functions by analysing
        the structure of the corresponding non-interacting Green's function.
        The results are stored in the
        :class:`SumkDFT.block_structure <dft.block_structure.BlockStructure>`
        attribute.

        Due to the implementation and numerics, the maximum difference between
        two matrix elements that are detected as equal can be a bit higher
        (e.g. a factor of two) than the actual threshold.

        Parameters
        ----------
        G : list of BlockGf of GfImFreq or GfImTime
            the non-interacting Green's function for each inequivalent correlated shell
        threshold : real, optional
                    If the difference between matrix elements is below threshold,
                    they are considered to be equal.
        include_shells : list of integers, optional
                         List of correlated shells to be analysed.
                         If include_shells is not provided all correlated shells will be analysed.
        """

        # initialize
        self.deg_shells = [[] for ish in range(self.n_inequiv_shells)]

        # helper function
        def null(A, eps=1e-15):
            """ Calculate the null-space of matrix A """
            u, s, vh = np.linalg.svd(A)
            null_mask = (s <= eps)
            null_space = compress(null_mask, vh, axis=0)
            return null_space.conjugate().transpose()

        gf = self._get_hermitian_quantity_from_gf(G)

        if include_shells is None:
            # include all shells
            include_shells = list(range(self.n_inequiv_shells))

        # We consider two blocks equal, if their Green's functions obey
        # maybe_conjugate1( v1^dagger G1 v1 ) = maybe_conjugate2( v2^dagger G2 v2 )
        # where maybe_conjugate is a function that conjugates the Green's
        # function if the flag 'conjugate' is set and the v are unitary
        # matrices
        #
        # for each pair of blocks, we check whether there is a transformation
        # maybe_conjugate( T G1 T^dagger ) = G2
        # where our goal is to find T
        # we just try whether there is such a T with and without conjugation
        for ish in include_shells:
            for block1 in self.gf_struct_solver[ish].keys():
                for block2 in self.gf_struct_solver[ish].keys():
                    if block1==block2: continue

                    # check if the blocks are already present in the deg_shells
                    ind1 = -1
                    ind2 = -2
                    for n, ind in enumerate(self.deg_shells[ish]):
                        if block1 in ind:
                            ind1 = n
                            v1 = ind[block1]
                        if block2 in ind:
                            ind2 = n
                            v2 = ind[block2]

                    # if both are already present, go to the next pair of blocks
                    if ind1 >= 0 and ind2 >= 0:
                        continue

                    gf1 = gf[ish][block1]
                    gf2 = gf[ish][block2]

                    # the two blocks have to have the same shape
                    if gf1.target_shape != gf2.target_shape:
                        continue

                    # Instead of directly comparing the two blocks, we
                    # compare its eigenvalues. As G(tau) is Hermitian,
                    # they are real and the eigenvector matrix is unitary.
                    # Thus, if the eigenvalues are equal we can transform
                    # one block to make it equal to the other (at least
                    # for tau=0).

                    e1 = np.linalg.eigvalsh(gf1.data[0])
                    e2 = np.linalg.eigvalsh(gf2.data[0])
                    if np.any(abs(e1-e2) > threshold): continue

                    for conjugate in [False,True]:
                        if conjugate:
                            gf2 = gf2.conjugate()

                        # we want T gf1 T^dagger = gf2
                        # while for a given tau, T could be calculated
                        # by diagonalizing gf1 and gf2, this does not
                        # work for all taus simultaneously because of
                        # numerical imprecisions

                        # rather, we rewrite the equation to
                        # T gf1 = gf2 T
                        # which is the Sylvester equation.
                        # For that equation, one can use the Kronecker
                        # product to get a linear problem, which consists
                        # of finding the null space of M vec T = 0.

                        M = np.kron(np.eye(*gf1.target_shape),gf2.data[0])-np.kron(gf1.data[0].transpose(),np.eye(*gf1.target_shape))
                        N = null(M, threshold)

                        # now we get the intersection of the null spaces
                        # of all values of tau
                        for i in range(1,len(gf1.data)):
                            M = np.kron(np.eye(*gf1.target_shape),gf2.data[i])-np.kron(gf1.data[i].transpose(),np.eye(*gf1.target_shape))
                            # transform M into current null space
                            M = np.dot(M, N)
                            N = np.dot(N, null(M, threshold))
                            if np.size(N) == 0:
                                break

                        # no intersection of the null spaces -> no symmetry
                        if np.size(N) == 0: continue

                        # reshape N: it then has the indices matrix, matrix, number of basis vectors of the null space
                        N = N.reshape(gf1.target_shape[0], gf1.target_shape[1], -1).transpose([1, 0, 2])

                        """
                        any matrix in the null space can now be constructed as
                        M = 0
                        for i in range(N.shape[-1]):
                            M += y[i]*N[:,:,i]
                        with coefficients (complex numbers) y[i].

                        We want to get a set of coefficients y so that M is unitary.
                        Unitary means M M^dagger = 1.
                        Thus,
                            sum  y[i] N[:,:,i] y[j].conjugate() N[:,:,j].conjugate().transpose() = eye.
                        The object N[:,:,i] N[:,:,j] is a four-index object which we call Z.
                        """
                        Z = np.einsum('aci,bcj->abij', N, N.conjugate())

                        """
                        function chi2
                        This function takes a real parameter vector y and reinterprets it as complex.
                        Then, it calculates the chi2 of
                            sum  y[i] N[:,:,i] y[j].conjugate() N[:,:,j].conjugate().transpose() - eye.
                        """
                        def chi2(y):
                            # reinterpret y as complex number
                            y = y.view(complex)
                            ret = 0.0
                            for a in range(Z.shape[0]):
                                for b in range(Z.shape[1]):
                                    ret += np.abs(np.dot(y, np.dot(Z[a, b], y.conjugate()))
                                                  - (1.0 if a == b else 0.0))**2
                            return ret

                        # use the minimization routine from scipy
                        res = minimize(chi2, np.ones(2 * N.shape[-1]))

                        # if the minimization fails, there is probably no symmetry
                        if not res.success: continue
                        # check if the minimization returned zero within the tolerance
                        if res.fun > threshold: continue

                        # reinterpret the solution as a complex number
                        y = res.x.view(complex)

                        # reconstruct the T matrix
                        T = np.zeros(N.shape[:-1], dtype=complex)
                        for i in range(len(y)):
                            T += N[:, :, i] * y[i]

                        # transform gf1 using T
                        G_transformed = gf1.copy()
                        G_transformed.from_L_G_R(T, gf1, T.conjugate().transpose())

                        # it does not make sense to check the tails for an
                        # absolute error because it will usually not hold;
                        # we could just check the relative error
                        # (here, we ignore it, reasoning that if the data
                        # is the same, the tails have to coincide as well)
                        try:
                            assert_arrays_are_close(G_transformed.data, gf2.data, threshold)
                        except (RuntimeError, AssertionError):
                            # the symmetry does not hold
                            continue

                        # Now that we have found a valid T, we have to
                        # rewrite it to match the convention that
                        # C1(v1^dagger G1 v1) = C2(v2^dagger G2 v2),
                        # where C conjugates if the flag is True

                        # For each group of degenerate shells, the list
                        # SK.deg_shells[ish] contains a dict. The keys
                        # of the dict are the block names, the values
                        # are tuples. The first entry of the tuple is
                        # the transformation matrix v, the second entry
                        # is the conjugation flag

                        # the second block is already present
                        # set v1 and C1 so that they are compatible with
                        # C(T gf1 T^dagger) = gf2
                        # and with
                        # C1(v1^dagger G1 v1) = C2(v2^dagger G2 v2)
                        if (ind1 < 0) and (ind2 >= 0):
                            if conjugate:
                                self.deg_shells[ish][ind2][block1] = np.dot(T.conjugate().transpose(), v2[0].conjugate()), not v2[1]
                            else:
                                self.deg_shells[ish][ind2][block1] = np.dot(T.conjugate().transpose(), v2[0]), v2[1]
                        # the first block is already present
                        # set v2 and C2 so that they are compatible with
                        # C(T gf1 T^dagger) = gf2
                        # and with
                        # C1(v1^dagger G1 v1) = C2(v2^dagger G2 v2)
                        elif (ind1 >= 0) and (ind2 < 0):
                            if conjugate:
                                self.deg_shells[ish][ind1][block2] = np.dot(T.conjugate(), v1[0].conjugate()), not v1[1]
                            else:
                                self.deg_shells[ish][ind1][block2] = np.dot(T, v1[0]), v1[1]
                        # the blocks are not already present
                        # we arbitrarily choose v1=eye and C1=False and
                        # set v2 and C2 so that they are compatible with
                        # C(T gf1 T^dagger) = gf2
                        # and with
                        # C1(v1^dagger G1 v1) = C2(v2^dagger G2 v2)
                        elif (ind1 < 0) and (ind2 < 0):
                            d = dict()
                            d[block1] = np.eye(*gf1.target_shape), False
                            if conjugate:
                                d[block2] = T.conjugate(), True
                            else:
                                d[block2] = T, False
                            self.deg_shells[ish].append(d)

                        # a block was found, break out of the loop
                        break

    def calculate_diagonalization_matrix(self, prop_to_be_diagonal='eal', calc_in_solver_blocks=True, write_to_blockstructure = True, shells=None):
        """
        Calculates the diagonalisation matrix, and (optionally) stores it in the BlockStructure.

        Parameters
        ----------
        prop_to_be_diagonal : string, optional
                              Defines the property to be diagonalized.

                              - 'eal' : local hamiltonian (i.e. crystal field)
                              - 'dm' : local density matrix

        calc_in_solver_blocks : bool, optional
                                Whether the property shall be diagonalized in the
                                full sumk structure, or just in the solver structure.

        write_to_blockstructure : bool, optional
                                Whether the diagonalization matrix shall be written to
                                the BlockStructure directly.
        shells : list of int, optional
              Indices of correlated shells to be diagonalized.
              None: all shells

        Returns
        -------
        trafo : dict
               The transformation matrix for each spin-block in the correlated shell
        """

        if self.block_structure.transformation:
            mpi.report(
                    "calculate_diagonalization_matrix: requires block_structure.transformation = None.")
            return 0

        # Use all shells
        if shells is None:
            shells = range(self.n_corr_shells)
        elif max(shells) >= self.n_corr_shells: # Check if the shell indices are present
            mpi.report("calculate_diagonalization_matrix: shells not correct.")
            return 0

        if prop_to_be_diagonal == 'eal':
            prop = [self.eff_atomic_levels()[self.corr_to_inequiv[ish]]
                    for ish in range(self.n_corr_shells)]
        elif prop_to_be_diagonal == 'dm':
            prop = self.density_matrix(method='using_point_integration')
        else:
            mpi.report(
                "calculate_diagonalization_matrix: Choices for prop_to_be_diagonal are 'eal' or 'dm'.")
            return 0

        trans = [{block: np.eye(block_dim) for block, block_dim in gfs} for gfs in self.gf_struct_sumk]

        for ish in shells:
            trafo = {}
            # Transform to solver basis if desired, blocks of prop change in this step!
            if calc_in_solver_blocks:
                prop[ish] = self.block_structure.convert_matrix(prop[ish], space_from='sumk', space_to='solver')
            # Get diagonalisation matrix, if not already diagonal
            for name in prop[ish]:
                if np.sum(abs(prop[ish][name]-np.diag(np.diagonal(prop[ish][name])))) > 1e-13:
                    trafo[name] = np.linalg.eigh(prop[ish][name])[1].conj().T
                else:
                    trafo[name] = np.identity(np.shape(prop[ish][name])[0])
            # Transform back to sumk if necessay, blocks change in this step!
            if calc_in_solver_blocks:
                trafo = self.block_structure.convert_matrix(trafo, space_from='solver', space_to='sumk')
            trans[ish] = trafo

        # Write to block_structure object

        if write_to_blockstructure:
            self.block_structure.transformation = trans

        return trans


    def density_matrix(self, method='using_gf'):
        """Calculate density matrices in one of two ways.

        Parameters
        ----------
        method : string, optional

                 - if 'using_gf': First get lattice gf (g_loc is not set up), then density matrix.
                                  It is useful for Hubbard I, and very quick.
                                  No assumption on the hopping structure is made (ie diagonal or not).
                 - if 'using_point_integration': Only works for diagonal hopping matrix (true in wien2k).

        beta : float, optional
               Inverse temperature.

        Returns
        -------
        dens_mat : list of dicts
                   Density matrix for each spin in each correlated shell.
        """
        dens_mat = [{} for icrsh in range(self.n_corr_shells)]
        for icrsh in range(self.n_corr_shells):
            for sp in self.spin_block_names[self.corr_shells[icrsh]['SO']]:
                dens_mat[icrsh][sp] = np.zeros(
                    [self.corr_shells[icrsh]['dim'], self.corr_shells[icrsh]['dim']], complex)

        ikarray = np.array(list(range(self.n_k)))
        for ik in mpi.slice_array(ikarray):

            if method == "using_gf":

                G_latt_iw = self.lattice_gf(ik=ik, mu=self.chemical_potential)
                G_latt_iw *= self.bz_weights[ik]
                dm = G_latt_iw.density()
                MMat = [dm[sp] for sp in self.spin_block_names[self.SO]]

            elif method == "using_point_integration":

                ntoi = self.spin_names_to_ind[self.SO]
                spn = self.spin_block_names[self.SO]
                dims = {sp:self.n_orbitals[ik, ntoi[sp]] for sp in spn}
                MMat = [np.zeros([dims[sp], dims[sp]], complex) for sp in spn]

                for isp, sp in enumerate(spn):
                    ind = ntoi[sp]
                    for inu in range(self.n_orbitals[ik, ind]):
                        # only works for diagonal hopping matrix (true in
                        # wien2k)
                        if (self.hopping[ik, ind, inu, inu] - self.h_field * (1 - 2 * isp)) < 0.0:
                            MMat[isp][inu, inu] = 1.0
                        else:
                            MMat[isp][inu, inu] = 0.0

            else:
                raise ValueError("density_matrix: the method '%s' is not supported." % method)

            for icrsh in range(self.n_corr_shells):
                for isp, sp in enumerate(self.spin_block_names[self.corr_shells[icrsh]['SO']]):
                    ind = self.spin_names_to_ind[
                        self.corr_shells[icrsh]['SO']][sp]
                    dim = self.corr_shells[icrsh]['dim']
                    n_orb = self.n_orbitals[ik, ind]
                    projmat = self.proj_mat[ik, ind, icrsh, 0:dim, 0:n_orb]
                    if method == "using_gf":
                        dens_mat[icrsh][sp] += np.dot(np.dot(projmat, MMat[isp]),
                                                         projmat.transpose().conjugate())
                    elif method == "using_point_integration":
                        dens_mat[icrsh][sp] += self.bz_weights[ik] * np.dot(np.dot(projmat, MMat[isp]),
                                                                               projmat.transpose().conjugate())

        # get data from nodes:
        for icrsh in range(self.n_corr_shells):
            for sp in dens_mat[icrsh]:
                dens_mat[icrsh][sp] = mpi.all_reduce(
                    mpi.world, dens_mat[icrsh][sp], lambda x, y: x + y)
        mpi.barrier()

        if self.symm_op != 0:
            dens_mat = self.symmcorr.symmetrize(dens_mat)

        # Rotate to local coordinate system:
        if self.use_rotations:
            for icrsh in range(self.n_corr_shells):
                for sp in dens_mat[icrsh]:
                    if self.rot_mat_time_inv[icrsh] == 1:
                        dens_mat[icrsh][sp] = dens_mat[icrsh][sp].conjugate()
                    dens_mat[icrsh][sp] = np.dot(np.dot(self.rot_mat[icrsh].conjugate().transpose(), dens_mat[icrsh][sp]),
                                                    self.rot_mat[icrsh])

        return dens_mat

    # For simple dft input, get crystal field splittings.
    def eff_atomic_levels(self):
        r"""
        Calculates the effective local Hamiltonian required as an input for
        the Hubbard I Solver.
        The local Hamiltonian (effective atomic levels) is calculated by
        projecting the on-site Bloch Hamiltonian:

        .. math:: H^{loc}_{m m'} = \sum_{k} P_{m \nu}(k) H_{\nu\nu'}(k) P^{*}_{\nu' m'}(k),

        where

        .. math:: H_{\nu\nu'}(k) = [\epsilon_{\nu k} - h_{z} \sigma_{z}] \delta_{\nu\nu'}.

        Parameters
        ----------
        None

        Returns
        -------
        eff_atlevels : gf_struct_sumk like
                       Effective local Hamiltonian :math:`H^{loc}_{m m'}` for each
                       inequivalent correlated shell.

        """

        # define matrices for inequivalent shells:
        eff_atlevels = [{} for ish in range(self.n_inequiv_shells)]
        for ish in range(self.n_inequiv_shells):
            for sp in self.spin_block_names[self.corr_shells[self.inequiv_to_corr[ish]]['SO']]:
                eff_atlevels[ish][sp] = np.identity(
                    self.corr_shells[self.inequiv_to_corr[ish]]['dim'], complex)
                eff_atlevels[ish][sp] *= -self.chemical_potential
                eff_atlevels[ish][
                    sp] -= self.dc_imp[self.inequiv_to_corr[ish]][sp]

        # sum over k:
        if not hasattr(self, "Hsumk"):
            # calculate the sum over k. Does not depend on mu, so do it only
            # once:
            self.Hsumk = [{} for icrsh in range(self.n_corr_shells)]
            for icrsh in range(self.n_corr_shells):
                dim = self.corr_shells[icrsh]['dim']
                for sp in self.spin_block_names[self.corr_shells[icrsh]['SO']]:
                    self.Hsumk[icrsh][sp] = np.zeros(
                        [dim, dim], complex)
                for isp, sp in enumerate(self.spin_block_names[self.corr_shells[icrsh]['SO']]):
                    ind = self.spin_names_to_ind[
                        self.corr_shells[icrsh]['SO']][sp]
                    for ik in range(self.n_k):
                        n_orb = self.n_orbitals[ik, ind]
                        MMat = np.identity(n_orb, complex)
                        MMat = self.hopping[
                            ik, ind, 0:n_orb, 0:n_orb] - (1 - 2 * isp) * self.h_field * MMat
                        projmat = self.proj_mat[ik, ind, icrsh, 0:dim, 0:n_orb]
                        self.Hsumk[icrsh][sp] += self.bz_weights[ik] * np.dot(np.dot(projmat, MMat),
                                                                                 projmat.conjugate().transpose())
            # symmetrisation:
            if self.symm_op != 0:
                self.Hsumk = self.symmcorr.symmetrize(self.Hsumk)

            # Rotate to local coordinate system:
            if self.use_rotations:
                for icrsh in range(self.n_corr_shells):
                    for sp in self.Hsumk[icrsh]:
                        if self.rot_mat_time_inv[icrsh] == 1:
                            self.Hsumk[icrsh][sp] = self.Hsumk[
                                icrsh][sp].conjugate()
                        self.Hsumk[icrsh][sp] = np.dot(np.dot(self.rot_mat[icrsh].conjugate().transpose(), self.Hsumk[icrsh][sp]),
                                                          self.rot_mat[icrsh])

        # add to matrix:
        for ish in range(self.n_inequiv_shells):
            for sp in eff_atlevels[ish]:
                eff_atlevels[ish][
                    sp] += self.Hsumk[self.inequiv_to_corr[ish]][sp]

        return eff_atlevels

    def init_dc(self):
        r"""
        Initializes the double counting terms.

        Parameters
        ----------
        None

        """
        self.dc_imp = [{} for icrsh in range(self.n_corr_shells)]
        for icrsh in range(self.n_corr_shells):
            dim = self.corr_shells[icrsh]['dim']
            spn = self.spin_block_names[self.corr_shells[icrsh]['SO']]
            for sp in spn:
                self.dc_imp[icrsh][sp] = np.zeros([dim, dim], float)
        self.dc_energ = [0.0 for icrsh in range(self.n_corr_shells)]

    def set_dc(self, dc_imp, dc_energ):
        r"""
        Sets double counting corrections to given values.

        Parameters
        ----------
        dc_imp : gf_struct_sumk like
                 Double-counting self-energy term.
        dc_energ : list of floats
                   Double-counting energy corrections for each correlated shell.

        """

        self.dc_imp = dc_imp
        self.dc_energ = dc_energ

    def calc_dc(self, dens_mat, orb=0, U_interact=None, J_hund=None,
                use_dc_formula=0, use_dc_value=None, transform=True, legacy=True):
        r"""
        Calculate and set the double counting corrections.

        If 'use_dc_value' is provided the double-counting term is uniformly initialized
        with this constant and 'U_interact' and 'J_hund' are ignored.

        If 'use_dc_value' is None the correction is evaluated according to
        one of the following formulae:

        * use_dc_formula = 0: fully-localised limit (FLL)
        * use_dc_formula = 1: Held's formula, i.e. mean-field formula for the Kanamori
                              type of the interaction Hamiltonian
        * use_dc_formula = 2: around mean-field (AMF)

        Note that FLL and AMF formulae were derived assuming a full Slater-type interaction
        term and should be thus used accordingly. For the Kanamori-type interaction
        one should use formula 1.

        The double-counting self-energy term is stored in `self.dc_imp` and the energy
        correction in `self.dc_energ`.

        Parameters
        ----------
        dens_mat : gf_struct_solver like
            Density matrix for the specified correlated shell.
        orb : int, optional
            Index of an inequivalent shell.
        U_interact : float, optional
            Value of interaction parameter `U`.
        J_hund : float, optional
            Value of interaction parameter `J`.
        use_dc_formula : int, optional
            Type of double-counting correction (see description).
        use_dc_value : float, optional
            Value of the double-counting correction. If specified
            `U_interact`, `J_hund` and `use_dc_formula` are ignored.
        transform : bool
            whether or not to use the transformation in block_structure
            to transform the dc
        """

        for icrsh in range(self.n_corr_shells):

            # ish is the index of the inequivalent shell corresponding to icrsh
            ish = self.corr_to_inequiv[icrsh]
            if ish != orb:
                continue  # ignore this orbital
            # *(1+self.corr_shells[icrsh]['SO'])
            dim = self.corr_shells[icrsh]['dim']
            spn = self.spin_block_names[self.corr_shells[icrsh]['SO']]

            Ncr = {sp: 0.0 for sp in spn}
            for block, inner in self.gf_struct_solver[ish].items():
                bl = self.solver_to_sumk_block[ish][block]
                Ncr[bl] += dens_mat[block].real.trace()
            Ncrtot = sum(Ncr.values())
            for sp in spn:
                self.dc_imp[icrsh][sp] = np.identity(dim, float)
                if self.SP == 0:  # average the densities if there is no SP:
                    Ncr[sp] = Ncrtot / len(spn)
                # correction for SO: we have only one block in this case, but
                # in DC we need N/2
                elif self.SP == 1 and self.SO == 1:
                    Ncr[sp] = Ncrtot / 2.0

            # Uses "orbital" dimension with SO for double counting
            if self.SP == 1 and self.SO == 1:
                assert dim % 2 == 0
                dim //= 2
            
            if use_dc_value is None:
                if legacy:

                    if U_interact is None and J_hund is None:
                        raise ValueError("set_dc: either provide U_interact and J_hund or set use_dc_value to dc value.")

                    if use_dc_formula == 0:  # FLL

                        self.dc_energ[icrsh] = U_interact / \
                            2.0 * Ncrtot * (Ncrtot - 1.0)
                        for sp in spn:
                            Uav = U_interact * (Ncrtot - 0.5) - \
                                J_hund * (Ncr[sp] - 0.5)
                            self.dc_imp[icrsh][sp] *= Uav
                            self.dc_energ[icrsh] -= J_hund / \
                                len(spn) * (Ncr[sp]) * (Ncr[sp] - 1.0)
                            mpi.report(
                                "DC for shell %(icrsh)i and block %(sp)s = %(Uav)f" % locals())

                    elif use_dc_formula == 1:  # Held's formula, with U_interact the interorbital onsite interaction

                        self.dc_energ[icrsh] = (U_interact + (dim - 1) * (U_interact - 2.0 * J_hund) + (
                            dim - 1) * (U_interact - 3.0 * J_hund)) / (2 * dim - 1) / 2.0 * Ncrtot * (Ncrtot - 1.0)
                        for sp in spn:
                            Uav = (U_interact + (dim - 1) * (U_interact - 2.0 * J_hund) + (dim - 1)
                                   * (U_interact - 3.0 * J_hund)) / (2 * dim - 1) * (Ncrtot - 0.5)
                            self.dc_imp[icrsh][sp] *= Uav
                            mpi.report(
                                "DC for shell %(icrsh)i and block %(sp)s = %(Uav)f" % locals())

                    elif use_dc_formula == 2:  # AMF

                        self.dc_energ[icrsh] = 0.5 * U_interact * Ncrtot * Ncrtot
                        for sp in spn:
                            Uav = U_interact * \
                                (Ncrtot - Ncr[sp] / dim) - \
                                J_hund * (Ncr[sp] - Ncr[sp] / dim)
                            self.dc_imp[icrsh][sp] *= Uav
                            self.dc_energ[
                                icrsh] -= (U_interact + (dim - 1) * J_hund) / dim / len(spn) * Ncr[sp] * Ncr[sp]
                            mpi.report(
                                "DC for shell %(icrsh)i and block %(sp)s = %(Uav)f" % locals())

                    mpi.report("DC energy for shell %s = %s" %
                               (icrsh, self.dc_energ[icrsh]))
                
                else:
                    #For legacy compatibility 
                    match use_dc_formula:
                        case 0:
                            mpi.report(f"Detected {use_dc_formula=}, changing to sFLL")
                            use_dc_formula = "sFLL"
                        case 1:
                            mpi.report(f"Detected {use_dc_formula=}, changing to cHeld")
                            use_dc_formula = "cHeld"
                        case 2:
                            mpi.report(f"Detected {use_dc_formula=}, changing to sAMF")
                            use_dc_formula = "sAMF"
                    
                    for sp in spn:
                        DC_val, E_val = compute_DC_from_density(N_tot=Ncrtot,U=U_interact, J=J_hund, n_orbitals=dim, N_spin=Ncr[sp], method=use_dc_formula)
                        self.dc_imp[icrsh][sp] *= DC_val
                        self.dc_energ[icrsh] -= E_val




            else:  # use value provided for user to determine dc_energ and dc_imp

                self.dc_energ[icrsh] = use_dc_value * Ncrtot
                for sp in spn:
                    self.dc_imp[icrsh][sp] *= use_dc_value

                mpi.report(
                    "DC for shell %(icrsh)i = %(use_dc_value)f" % locals())
                mpi.report("DC energy = %s" % self.dc_energ[icrsh])
            if transform:
                for sp in spn:
                    T = self.block_structure.effective_transformation_sumk[icrsh][sp]
                    self.dc_imp[icrsh][sp] = np.dot(T.conjugate().transpose(),
                            np.dot(self.dc_imp[icrsh][sp], T))

    def add_dc(self):
        r"""
        Subtracts the double counting term from the impurity self energy.

        Parameters
        ----------

        Returns
        -------
        sigma_minus_dc : gf_struct_sumk like
                         Self-energy with a subtracted double-counting term.

        """

        # Be careful: Sigma_imp is already in the global coordinate system!!
        sigma_minus_dc = [s.copy() for s in self.Sigma_imp]
        for icrsh in range(self.n_corr_shells):
            for bname, gf in sigma_minus_dc[icrsh]:
                # Transform dc_imp to global coordinate system
                if self.use_rotations:
                    gf -= np.dot(self.rot_mat[icrsh], np.dot(self.dc_imp[icrsh][
                                       bname], self.rot_mat[icrsh].conjugate().transpose()))
                else:
                    gf -= self.dc_imp[icrsh][bname]

        return sigma_minus_dc

    def symm_deg_gf(self, gf_to_symm, ish=0):
        r"""
        Averages a GF over degenerate shells.

        Degenerate shells of an inequivalent correlated shell are defined by
        `self.deg_shells`. This function enforces corresponding degeneracies
        in the input GF.

        Parameters
        ----------
        gf_to_symm : gf_struct_solver like
                     Input and output GF (i.e., it gets overwritten)
        ish : int
              Index of an inequivalent shell. (default value 0)

        """

        # when reading block_structures written with older versions from
        # an h5 file, self.deg_shells might be None
        if self.deg_shells is None: return

        for degsh in self.deg_shells[ish]:
            # ss will hold the averaged orbitals in the basis where the
            # blocks are all equal
            # i.e. maybe_conjugate(v^dagger gf v)
            ss = None
            n_deg = len(degsh)
            for key in degsh:
                if ss is None:
                    ss = gf_to_symm[key].copy()
                    ss.zero()
                    helper = ss.copy()
                # get the transformation matrix
                if isinstance(degsh, dict):
                    v, C = degsh[key]
                else:
                    # for backward compatibility, allow degsh to be a list
                    v = np.eye(*ss.target_shape)
                    C = False
                # the helper is in the basis where the blocks are all equal
                helper.from_L_G_R(v.conjugate().transpose(), gf_to_symm[key], v)
                if C:
                    helper << helper.transpose()
                # average over all shells
                ss += helper / (1.0 * n_deg)
            # now put back the averaged gf to all shells
            for key in degsh:
                if isinstance(degsh, dict):
                    v, C = degsh[key]
                else:
                    # for backward compatibility, allow degsh to be a list
                    v = np.eye(*ss.target_shape)
                    C = False
                if C:
                    gf_to_symm[key].from_L_G_R(v, ss.transpose().copy(), v.conjugate().transpose())
                else:
                    gf_to_symm[key].from_L_G_R(v, ss, v.conjugate().transpose())

    def total_density(self, mu=None, with_Sigma=True, with_dc=True, broadening=None):
        r"""
        Calculates the total charge within the energy window for a given chemical potential.
        The chemical potential is either given by parameter `mu` or, if it is not specified,
        taken from `self.chemical_potential`.

        The total charge is calculated from the trace of the GF in the Bloch basis.
        By default, a full interacting GF is used. To use the non-interacting GF, set
        parameter `with_Sigma = False`.

        The number of bands within the energy windows generally depends on `k`. The trace is
        therefore calculated separately for each `k`-point.

        Since in general n_orbitals depends on k, the calculation is done in the following order:

        .. math:: n_{tot} = \sum_{k} n(k),

        with

        .. math:: n(k) = Tr G_{\nu\nu'}(k, i\omega_{n}).

        The calculation is done in the global coordinate system, if distinction is made between local/global.

        Parameters
        ----------
        mu : float, optional
             Input chemical potential. If not specified, `self.chemical_potential` is used instead.
        with_Sigma : boolean, optional
             If `True` the full interacing GF is evaluated, otherwise the self-energy is not
             included and the charge would correspond to a non-interacting system.
        with_dc : boolean, optional
             Whether or not to subtract the double-counting term from the self-energy.
        broadening : float, optional
                     Imaginary shift for the axis along which the real-axis GF is calculated.
                     If not provided, broadening will be set to double of the distance between mesh points in 'mesh'.
                     Only relevant for real-frequency GF.

        Returns
        -------
        dens : float
               Total charge :math:`n_{tot}`.

        """

        if mu is None:
            mu = self.chemical_potential
        dens = 0.0
        ikarray = np.array(list(range(self.n_k)))
        for ik in mpi.slice_array(ikarray):
            G_latt = self.lattice_gf(
                ik=ik, mu=mu, with_Sigma=with_Sigma, with_dc=with_dc, broadening=broadening)
            dens += self.bz_weights[ik] * G_latt.total_density()
        # collect data from mpi:
        dens = mpi.all_reduce(mpi.world, dens, lambda x, y: x + y)
        mpi.barrier()

        if abs(dens.imag) > 1e-20:
            mpi.report("Warning: Imaginary part in density will be ignored ({})".format(str(abs(dens.imag))))
        return dens.real

    def set_mu(self, mu):
        r"""
        Sets a new chemical potential.

        Parameters
        ----------
        mu : float
             New value of the chemical potential.

        """
        self.chemical_potential = mu

    def calc_mu(self, precision=0.01, broadening=None, delta=0.5, max_loops=100, method="dichotomy"):
        r"""
        Searches for the chemical potential that gives the DFT total charge.
        A simple bisection method is used.

        Parameters
        ----------
        precision : float, optional
                    A desired precision of the resulting total charge.
        broadening : float, optional
                     Imaginary shift for the axis along which the real-axis GF is calculated.
                     If not provided, broadening will be set to double of the distance between mesh points in 'mesh'.
                     Only relevant for real-frequency GF.
        max_loops : int, optional
                    Number of dichotomy loops maximally performed.
        
        max_loops : string, optional
                    Type of optimization used:
                        * dichotomy: usual bisection algorithm from the TRIQS library
                        * newton: newton method, faster convergence but more unstable 
                        * brent: finds bounds and proceeds with hyperbolic brent method, a compromise between speed and ensuring convergence

        Returns
        -------
        mu : float
             Value of the chemical potential giving the DFT total charge
             within specified precision.

        """
        def find_bounds(function, x_init, delta_x, max_loops=1000):
            mpi.report("Finding bounds on chemical potential" )
            x= x_init
            # First find the bounds
            y1 = function(x)
            eps = np.sign(y1)
            x1= x;
            x2= x1;y2 = y1
            
            nbre_loop=0
            #abort the loop after maxiter is reached or when y1 and y2 have different sign 
            while (nbre_loop<= max_loops) and (y2*y1)>0:
                nbre_loop +=1
                x1=x2
                y1=y2

                x2 -=  eps*delta_x
                y2 = function(x2)

            if nbre_loop > (max_loops):
                raise ValueError("The bounds could not be found")

            # Make sure that x2 > x1
            if x1 > x2:
                x1,x2 = x2,x1
                y1,y2 = y2,y1


            mpi.report(f"mu_interval: [  {x1:.4f}  ; {x2:.4f} ]")
            mpi.report(f"delta to target density interval: [ {y1:.4f} ; {y2:.4f} ]")
            return x1, x2



        # previous implementation
        def F_bisection(mu): return self.total_density(mu=mu, broadening=broadening).real
        density = self.density_required - self.charge_below
        #using scipy.optimize
        def F_optimize(mu):

            mpi.report("Trying out mu = {}".format(str(mu)))
            calc_dens = self.total_density(mu=mu, broadening=broadening).real - density
            mpi.report("Delta to target density = {}".format(str(calc_dens)))
            return calc_dens
        
        #check for lowercase matching for the method variable
        match method.lower():

            case "dichotomy":
                mpi.report("\nSUMK calc_mu: Using dichtomy adjustment to find chemical potential\n")
                self.chemical_potential = dichotomy.dichotomy(function=F_bisection,
                                                              x_init=self.chemical_potential, y_value=density,
                                                              precision_on_y=precision, delta_x=delta, max_loops=max_loops,
                                                              x_name="Chemical Potential", y_name="Total Density",
                                                              verbosity=3)[0]
            case 'newton':
                mpi.report("\nSUMK calc_mu: Using Newton method to find chemical potential\n")
                self.chemical_potential =                newton(func=F_optimize,
                                                              x0=self.chemical_potential,
                                                              tol=precision, maxiter=max_loops,
                                                              )

            case 'brent':
                mpi.report("\nSUMK calc_mu: Using Brent method to find chemical potential")
                mpi.report("SUMK calc_mu: Finding bounds \n")
                
                mu_guess_0, mu_guess_1  =  find_bounds(function=F_optimize,
                                                      x_init=self.chemical_potential,
                                                      delta_x=delta, max_loops=max_loops,
                                                      )
                mu_guess_1 += 0.01 #scrambles higher lying interval to avoid getting stuck
                mpi.report("\nSUMK calc_mu: Searching root with Brent method\n")
                self.chemical_potential = brenth(f=F_optimize,
                                                   a=mu_guess_0,
                                                   b=mu_guess_1,
                                                   xtol=precision, maxiter=max_loops,
                                                   )
            
            case _:
                raise ValueError(
                        f"SUMK calc_mu: The method selected: {method}, is not implemented\n",
                        """
                        Please check for typos or select one of the following:
                            * dichotomy: usual bisection algorithm from the TRIQS library
                            * newton: newton method, fastest convergence but more unstable 
                            * brent: finds bounds and proceeds with hyperbolic brent method, a compromise between speed and ensuring convergence
                        """
                        )

        return self.chemical_potential

    def calc_density_correction(self, filename=None, dm_type=None, spinave=False, kpts_to_write=None):
        r"""
        Calculates the charge density correction and stores it into a file.

        The charge density correction is needed for charge-self-consistent DFT+DMFT calculations.
        It represents a density matrix of the interacting system defined in Bloch basis
        and it is calculated from the sum over Matsubara frequecies of the full GF,

        ..math:: N_{\nu\nu'}(k) = \sum_{i\omega_{n}} G_{\nu\nu'}(k, i\omega_{n})

        The density matrix for every `k`-point is stored into a file.

        Parameters
        ----------
        filename : string
                   Name of the file to store the charge density correction.
        dm_type : string
                   DFT code to write the density correction for. Options:
                   'vasp', 'wien2k', 'elk' or 'qe'. Needs to be set for 'qe'
        spinave : logical
                   Elk specific and for magnetic calculations in DMFT only. 
                   It averages the spin to keep the DFT part non-magnetic.            
        kpts_to_write : iterable of int
                   Indices of k points that are written to file. If None (default),
                   all k points are written. Only implemented for dm_type 'vasp'

        Returns
        -------
        (deltaN, dens) : tuple
                         Returns a tuple containing the density matrix `deltaN` and
                         the corresponing total charge `dens`.

        """
        #automatically set dm_type if required
        if dm_type==None:
            dm_type = self.dft_code

        assert dm_type in ('vasp', 'wien2k','elk', 'qe'), "'dm_type' must be either 'vasp', 'wienk', 'elk' or 'qe'"
        #default file names
        if filename is None:
            if dm_type == 'wien2k':
                filename = 'dens_mat.dat'
            elif dm_type == 'vasp':
                filename = 'GAMMA'
            elif dm_type == 'elk':
                filename = 'DMATDMFT.OUT'
            elif dm_type == 'qe':
                filename = self.hdf_file


        assert isinstance(filename, str), ("calc_density_correction: "
                                              "filename has to be a string!")

        assert kpts_to_write is None or dm_type == 'vasp', ('Selecting k-points only'
                                                            +'implemented for vasp')

        ntoi = self.spin_names_to_ind[self.SO]
        spn = self.spin_block_names[self.SO]
        dens = {sp: 0.0 for sp in spn}
        band_en_correction = 0.0

# Fetch Fermi weights and energy window band indices
        if dm_type in ['vasp','qe']:
            fermi_weights = 0
            band_window = 0
            if mpi.is_master_node():
                with HDFArchive(self.hdf_file,'r') as ar:
                    fermi_weights = ar['dft_misc_input']['dft_fermi_weights']
                    band_window = ar['dft_misc_input']['band_window']
            fermi_weights = mpi.bcast(fermi_weights)
            band_window = mpi.bcast(band_window)

# Convert Fermi weights to a density matrix
            dens_mat_dft = {}
            for sp in spn:
                dens_mat_dft[sp] = [fermi_weights[ik, ntoi[sp], :].astype(complex) for ik in range(self.n_k)]


        # Set up deltaN:
        deltaN = {}
        for sp in spn:
            deltaN[sp] = [np.zeros([self.n_orbitals[ik, ntoi[sp]], self.n_orbitals[
                                      ik, ntoi[sp]]], complex) for ik in range(self.n_k)]

        ikarray = np.arange(self.n_k)
        for ik in mpi.slice_array(ikarray):
            G_latt_iw = self.lattice_gf(
                ik=ik, mu=self.chemical_potential)
            if dm_type == 'vasp' and self.proj_or_hk == 'hk':
                # rotate the Green function into the DFT band basis
                for bname, gf in G_latt_iw:
                    G_latt_rot_iw = gf.copy()
                    G_latt_rot_iw << self.upfold(
                            ik, 0, bname, G_latt_iw[bname], gf,shells='csc')

                    G_latt_iw[bname] = G_latt_rot_iw.copy()

            for bname, gf in G_latt_iw:
                deltaN[bname][ik] = G_latt_iw[bname].density()

                dens[bname] += self.bz_weights[ik] * G_latt_iw[bname].total_density()
                if dm_type in ['vasp','qe']:
# In 'vasp'-mode subtract the DFT density matrix
                    nb = self.n_orbitals[ik, ntoi[bname]]
                    diag_inds = np.diag_indices(nb)
                    deltaN[bname][ik][diag_inds] -= dens_mat_dft[bname][ik][:nb]

                    if self.charge_mixing and self.deltaNOld is not None:
                        G2 = np.sum(self.kpts_cart[ik,:]**2)
                        # Kerker mixing
                        mix_fac = self.charge_mixing_alpha * G2 / (G2 + self.charge_mixing_gamma**2)
                        deltaN[bname][ik][diag_inds] = (1.0 - mix_fac) * self.deltaNOld[bname][ik][diag_inds] + mix_fac * deltaN[bname][ik][diag_inds]
                    dens[bname] -= self.bz_weights[ik] * dens_mat_dft[bname][ik].sum().real
                    isp = ntoi[bname]
                    b1, b2 = band_window[isp][ik, :2]
                    nb = b2 - b1 + 1
                    assert nb == self.n_orbitals[ik, ntoi[bname]], "Number of bands is inconsistent at ik = %s"%(ik)
                    band_en_correction += np.dot(deltaN[bname][ik], self.hopping[ik, isp, :nb, :nb]).trace().real * self.bz_weights[ik]

        # mpi reduce:
        for bname in deltaN:
            for ik in range(self.n_k):
                deltaN[bname][ik] = mpi.all_reduce(
                    mpi.world, deltaN[bname][ik], lambda x, y: x + y)
            dens[bname] = mpi.all_reduce(
                mpi.world, dens[bname], lambda x, y: x + y)
        self.deltaNOld = copy.copy(deltaN)
        mpi.barrier()



        band_en_correction = mpi.all_reduce(mpi.world, band_en_correction, lambda x,y : x+y)

        # now save to file:
        if dm_type == 'wien2k':
            if mpi.is_master_node():
                if self.SP == 0:
                    f = open(filename, 'w')
                else:
                    f = open(filename + 'up', 'w')
                    f1 = open(filename + 'dn', 'w')
                # write chemical potential (in Rydberg):
                f.write("%.14f\n" % (self.chemical_potential / self.energy_unit))
                if self.SP != 0:
                    f1.write("%.14f\n" %
                             (self.chemical_potential / self.energy_unit))
                # write beta in rydberg-1
                f.write("%.14f\n" % (G_latt_iw.mesh.beta * self.energy_unit))
                if self.SP != 0:
                    f1.write("%.14f\n" % (G_latt_iw.mesh.beta * self.energy_unit))

                if self.SP == 0:  # no spin-polarization

                    for ik in range(self.n_k):
                        f.write("%s\n" % self.n_orbitals[ik, 0])
                        for inu in range(self.n_orbitals[ik, 0]):
                            for imu in range(self.n_orbitals[ik, 0]):
                                valre = (deltaN['up'][ik][
                                         inu, imu].real + deltaN['down'][ik][inu, imu].real) / 2.0
                                valim = (deltaN['up'][ik][
                                         inu, imu].imag + deltaN['down'][ik][inu, imu].imag) / 2.0
                                f.write("%.14f  %.14f " % (valre, valim))
                            f.write("\n")
                        f.write("\n")
                    f.close()

                elif self.SP == 1:  # with spin-polarization

                    # dict of filename: (spin index, block_name)
                    if self.SO == 0:
                        to_write = {f: (0, 'up'), f1: (1, 'down')}
                    if self.SO == 1:
                        to_write = {f: (0, 'ud'), f1: (0, 'ud')}
                    for fout in to_write.keys():
                        isp, sp = to_write[fout]
                        for ik in range(self.n_k):
                            fout.write("%s\n" % self.n_orbitals[ik, isp])
                            for inu in range(self.n_orbitals[ik, isp]):
                                for imu in range(self.n_orbitals[ik, isp]):
                                    fout.write("%.14f  %.14f " % (deltaN[sp][ik][
                                               inu, imu].real, deltaN[sp][ik][inu, imu].imag))
                                fout.write("\n")
                            fout.write("\n")
                        fout.close()
        elif dm_type == 'vasp':
            if kpts_to_write is None:
                kpts_to_write = np.arange(self.n_k)
            else:
                assert np.min(kpts_to_write) >= 0 and np.max(kpts_to_write) < self.n_k

            assert self.SP == 0, "Spin-polarized density matrix is not implemented"

            if mpi.is_master_node():
                with open(filename, 'w') as f:
                    f.write(" %i  -1  ! Number of k-points, default number of bands\n"%len(kpts_to_write))
                    for index, ik in enumerate(kpts_to_write):
                        ib1 = band_window[0][ik, 0]
                        ib2 = band_window[0][ik, 1]
                        f.write(" %i  %i  %i\n"%(index + 1, ib1, ib2))
                        for inu in range(self.n_orbitals[ik, 0]):
                            for imu in range(self.n_orbitals[ik, 0]):
                                valre = (deltaN['up'][ik][inu, imu].real + deltaN['down'][ik][inu, imu].real) / 2.0
                                valim = (deltaN['up'][ik][inu, imu].imag + deltaN['down'][ik][inu, imu].imag) / 2.0
                                f.write(" %.14f  %.14f"%(valre, valim))
                            f.write("\n")

        elif dm_type == 'elk':
        # output each k-point density matrix for Elk
            if mpi.is_master_node():
        # read in misc data from .h5 file
                things_to_read = ['band_window','vkl','nstsv']
                self.subgroup_present, self.value_read = self.read_input_from_hdf(
                             subgrp=self.misc_data, things_to_read=things_to_read)
        # open file
                with open(filename, 'w') as f:
        # determine the number of spin blocks
                    n_spin_blocks = self.SP + 1 - self.SO
                    nbmax = np.max(self.n_orbitals)
        # output beta and mu in Hartrees
                    beta = G_latt_iw.mesh.beta * self.energy_unit
                    mu = self.chemical_potential/self.energy_unit
        # ouput n_k, nspin and max orbitals - a check
                    f.write(" %d  %d  %d  %.14f %.14f ! nkpt, nspin, nstmax, beta, mu\n"%(self.n_k, n_spin_blocks, nbmax, beta, mu))
                    for ik in range(self.n_k):
                      for ispn in range(n_spin_blocks):
                        #Determine the SO density matrix band indices from the spinor band indices
                        if(self.SO==1):
                          band0=[self.band_window[0][ik, 0],self.band_window[1][ik, 0]]
                          band1=[self.band_window[0][ik, 1],self.band_window[1][ik, 1]]
                          ib1=int(min(band0))
                          ib2=int(max(band1))
                        else:
                        #Determine the density matrix band indices from the spinor band indices
                          ib1 = self.band_window[ispn][ik, 0]
                          ib2 = self.band_window[ispn][ik, 1]
                        f.write(" %d  %d  %d  %d ! ik, ispn, minist, maxist\n"%(ik + 1, ispn + 1, ib1, ib2))
                        for inu in range(self.n_orbitals[ik, ispn]):
                            for imu in range(self.n_orbitals[ik, ispn]):
                                #output non-magnetic or spin-averaged density matrix
                                if((self.SP==0) or (spinave)):
                                  valre = (deltaN['up'][ik][inu, imu].real + deltaN['down'][ik][inu, imu].real) / 2.0
                                  valim = (deltaN['up'][ik][inu, imu].imag + deltaN['down'][ik][inu, imu].imag) / 2.0
                                else:
                                  valre = deltaN[spn[ispn]][ik][inu, imu].real
                                  valim = deltaN[spn[ispn]][ik][inu, imu].imag
                                f.write(" %.14f  %.14f"%(valre, valim))
                            f.write("\n")

        elif dm_type == 'qe':
            assert self.SP == 0, "Spin-polarized density matrix is not implemented"

            subgrp = 'dft_update'
            delta_N = np.zeros([self.n_k, max(self.n_orbitals[:,0]), max(self.n_orbitals[:,0])], dtype=complex)
            mpi.report(" %i  -1  ! Number of k-points, default number of bands\n"%(self.n_k))
            for ik in range(self.n_k):
                ib1 = band_window[0][ik, 0]
                ib2 = band_window[0][ik, 1]
                for inu in range(self.n_orbitals[ik, 0]):
                    for imu in range(self.n_orbitals[ik, 0]):
                        valre = (deltaN['up'][ik][inu, imu].real + deltaN['down'][ik][inu, imu].real) / 2.0
                        valim = (deltaN['up'][ik][inu, imu].imag + deltaN['down'][ik][inu, imu].imag) / 2.0
                        # write into delta_N
                        delta_N[ik, inu, imu] = valre + 1j*valim
            if mpi.is_master_node():
                with HDFArchive(self.hdf_file, 'a') as ar:
                    if not subgrp in ar:
                        ar.create_group(subgrp)
                    things_to_save = ['delta_N']
                    for it in things_to_save:
                        ar[subgrp][it] = locals()[it]


        else:
            raise NotImplementedError("Unknown density matrix type: '%s'"%(dm_type))

        res = deltaN, dens

        if dm_type in ['vasp', 'qe']:
            res += (band_en_correction,)

        return res

    def calculate_min_max_band_energies(self):
        hop = self.hopping
        diag_hop = np.zeros(hop.shape[:-1])
        hop_slice = mpi.slice_array(hop)
        diag_hop_slice = mpi.slice_array(diag_hop)
        diag_hop_slice[:] = np.linalg.eigvalsh(hop_slice)
        diag_hop = mpi.all_reduce(mpi.world, diag_hop, lambda x, y: x + y)
        min_band_energy = diag_hop.min().real
        max_band_energy = diag_hop.max().real
        self.min_band_energy = min_band_energy
        self.max_band_energy = max_band_energy
        return min_band_energy, max_band_energy

################
# FIXME LEAVE UNDOCUMENTED
################

    def calc_dc_for_density(self, orb, dc_init, dens_mat, density=None, precision=0.01):
        """Searches for DC in order to fulfill charge neutrality.
           If density is given, then DC is set such that the LOCAL charge of orbital
           orb coincides with the given density."""

        def F(dc):
            self.calc_dc(dens_mat=dens_mat, U_interact=0,
                         J_hund=0, orb=orb, use_dc_value=dc)
            if dens_req is None:
                return self.total_density(mu=mu)
            else:
                return self.extract_G_loc()[orb].total_density()

        if density is None:
            density = self.density_required - self.charge_below

        dc = dichotomy.dichotomy(function=F,
                                 x_init=dc_init, y_value=density,
                                 precision_on_y=precision, delta_x=0.5,
                                 max_loops=100, x_name="Double Counting", y_name="Total Density",
                                 verbosity=3)[0]

        return dc

    def check_projectors(self):
        """Calculated the density matrix from projectors (DM = P Pdagger) to check that it is correct and
           specifically that it matches DFT."""
        dens_mat = [np.zeros([self.corr_shells[icrsh]['dim'], self.corr_shells[icrsh]['dim']], complex)
                    for icrsh in range(self.n_corr_shells)]

        for ik in range(self.n_k):
            for icrsh in range(self.n_corr_shells):
                dim = self.corr_shells[icrsh]['dim']
                n_orb = self.n_orbitals[ik, 0]
                projmat = self.proj_mat[ik, 0, icrsh, 0:dim, 0:n_orb]
                dens_mat[icrsh][
                    :, :] += np.dot(projmat, projmat.transpose().conjugate()) * self.bz_weights[ik]

        if self.symm_op != 0:
            dens_mat = self.symmcorr.symmetrize(dens_mat)

        # Rotate to local coordinate system:
        if self.use_rotations:
            for icrsh in range(self.n_corr_shells):
                if self.rot_mat_time_inv[icrsh] == 1:
                    dens_mat[icrsh] = dens_mat[icrsh].conjugate()
                dens_mat[icrsh] = np.dot(np.dot(self.rot_mat[icrsh].conjugate().transpose(), dens_mat[icrsh]),
                                            self.rot_mat[icrsh])

        return dens_mat

    def sorts_of_atoms(self, shells):
        """
        Determine the number of inequivalent sorts.
        """
        sortlst = [shells[i]['sort'] for i in range(len(shells))]
        n_sorts = len(set(sortlst))
        return n_sorts

    def number_of_atoms(self, shells):
        """
        Determine the number of inequivalent atoms.
        """
        atomlst = [shells[i]['atom'] for i in range(len(shells))]
        n_atoms = len(set(atomlst))
        return n_atoms

    # The following methods are here to ensure backward-compatibility
    # after introducing the block_structure class
    def __get_gf_struct_sumk(self):
        return self.block_structure.gf_struct_sumk
    def __set_gf_struct_sumk(self,value):
        self.block_structure.gf_struct_sumk = value
    gf_struct_sumk = property(__get_gf_struct_sumk,__set_gf_struct_sumk)

    def __get_gf_struct_solver(self):
        return self.block_structure.gf_struct_solver
    def __set_gf_struct_solver(self,value):
        self.block_structure.gf_struct_solver = value
    gf_struct_solver = property(__get_gf_struct_solver,__set_gf_struct_solver)

    def __get_solver_to_sumk(self):
        return self.block_structure.solver_to_sumk
    def __set_solver_to_sumk(self,value):
        self.block_structure.solver_to_sumk = value
    solver_to_sumk = property(__get_solver_to_sumk,__set_solver_to_sumk)

    def __get_sumk_to_solver(self):
        return self.block_structure.sumk_to_solver
    def __set_sumk_to_solver(self,value):
        self.block_structure.sumk_to_solver = value
    sumk_to_solver = property(__get_sumk_to_solver,__set_sumk_to_solver)

    def __get_solver_to_sumk_block(self):
        return self.block_structure.solver_to_sumk_block
    def __set_solver_to_sumk_block(self,value):
        self.block_structure.solver_to_sumk_block = value
    solver_to_sumk_block = property(__get_solver_to_sumk_block,__set_solver_to_sumk_block)

    def __get_deg_shells(self):
        return self.block_structure.deg_shells
    def __set_deg_shells(self,value):
        self.block_structure.deg_shells = value
    deg_shells = property(__get_deg_shells,__set_deg_shells)

    @property
    def gf_struct_solver_list(self):
        return self.block_structure.gf_struct_solver_list

    @property
    def gf_struct_sumk_list(self):
        return self.block_structure.gf_struct_sumk_list

    @property
    def gf_struct_solver_dict(self):
        return self.block_structure.gf_struct_solver_dict

    @property
    def gf_struct_sumk_dict(self):
        return self.block_structure.gf_struct_sumk_dict

    def __get_corr_to_inequiv(self):
        return self.block_structure.corr_to_inequiv
    def __set_corr_to_inequiv(self, value):
        self.block_structure.corr_to_inequiv = value
    corr_to_inequiv = property(__get_corr_to_inequiv, __set_corr_to_inequiv)

    def __get_inequiv_to_corr(self):
        return self.block_structure.inequiv_to_corr
    def __set_inequiv_to_corr(self, value):
        self.block_structure.inequiv_to_corr = value
    inequiv_to_corr = property(__get_inequiv_to_corr, __set_inequiv_to_corr)
