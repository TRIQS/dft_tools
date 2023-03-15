##########################################################################

#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2011 by M. Aichhorn, L. Pourovskii, V. Vildosola
# Copyright (c) 2022-2023 Simons Foundation
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
# Authors: M. Aichhorn, S. Beck, A. Hampel, L. Pourovskii, V. Vildosola

##########################################################################
import sys
import numpy
from warnings import warn
from triqs.gf import *
import triqs.utility.mpi as mpi
from .symmetry import *
import scipy.constants as cst
import os.path

__all__ = ['transport_distribution', 'conductivity_and_seebeck', 'write_output_to_hdf',
           'init_spectroscopy', 'transport_function']

# ----------------- helper functions -----------------------


def read_transport_input_from_hdf(sum_k):
    r"""
    Reads the data for transport calculations from the hdf5 archive.
    Parameters
    ----------
    sum_k : sum_k object
            triqs SumkDFT object
    Returns
    -------
    sum_k : sum_k object
            triqs SumkDFT object
    """

    assert sum_k.dft_code in (
        'wien2k', 'elk', 'w90'), "read_transport_input_from_hdf() is only implemented for wien2k and elk inputs"

    if sum_k.dft_code in ('wien2k', 'elk'):
        thingstoread = ['band_window_optics', 'velocities_k']
    else:
        thingstoread = ['band_window_optics']

    sum_k.read_input_from_hdf(subgrp=sum_k.transp_data, things_to_read=thingstoread)

    if (sum_k.dft_code == "wien2k"):
        thingstoread = ['band_window', 'lattice_angles', 'lattice_constants',
                        'lattice_type', 'n_symmetries', 'rot_symmetries']
    elif (sum_k.dft_code == "elk"):
        thingstoread = ['band_window', 'n_symmetries', 'rot_symmetries',
                        'cell_vol']
    elif (sum_k.dft_code == 'w90'):
        thingstoread = ['band_window', 'n_symmetries', 'rot_symmetries']

    sum_k.read_input_from_hdf(subgrp=sum_k.misc_data, things_to_read=thingstoread)

    if (sum_k.dft_code == "wien2k"):
        sum_k.cell_vol = cellvolume(
            sum_k.lattice_type, sum_k.lattice_constants, sum_k.lattice_angles)[1]

    return sum_k


def write_output_to_hdf(sum_k, things_to_save, subgrp='user_data'):
    r"""
    Saves data from a list into the HDF file. Prints a warning if a requested data is not found in SumkDFT object.

    Parameters
    ----------
    hdf_file : hdf5 archive
               hd5 file
    things_to_save : list of strings
                     List of datasets to be saved into the hdf5 file.
    subgrp : string, optional
             Name of hdf5 file subgroup in which the data are to be stored.
    """

    if not (mpi.is_master_node()):
        return  # do nothing on nodes
    with HDFArchive(sum_k.hdf_file, 'a') as ar:
        if not subgrp in ar:
            ar.create_group(subgrp)
        for it, val in things_to_save.items():
            if it in ["gf_struct_sumk", "gf_struct_solver",
                      "solver_to_sumk", "sumk_to_solver", "solver_to_sumk_block"]:
                warn("It is not recommended to save '{}' individually. Save 'block_structure' instead.".format(it))
            ar[subgrp][it] = val


def cellvolume(lattice_type, lattice_constants, latticeangle):
    r"""
    Determines the conventional und primitive unit cell volumes.

    Parameters
    ----------
    lattice_type : string
        Lattice type according to the Wien2k convention (P, F, B, R, H, CXY, CYZ, CXZ).
    lattice_constants : list of double
        Lattice constants (a, b, c).
    lattice angles : list of double
        Lattice angles (:math:`\alpha, \beta, \gamma`).

    Returns
    -------
    vol_c : double
        Conventional unit cell volume.
    vol_p : double
        Primitive unit cell volume.
    """

    a = lattice_constants[0]
    b = lattice_constants[1]
    c = lattice_constants[2]
    c_al = numpy.cos(latticeangle[0])
    c_be = numpy.cos(latticeangle[1])
    c_ga = numpy.cos(latticeangle[2])
    vol_c = a * b * c * \
        numpy.sqrt(1 + 2 * c_al * c_be * c_ga -
                   c_al ** 2 - c_be ** 2 - c_ga ** 2)

    det = {"P": 1, "F": 4, "B": 2, "R": 3,
           "H": 1, "CXY": 2, "CYZ": 2, "CXZ": 2}
    vol_p = vol_c / det[lattice_type]

    return vol_c, vol_p


def fermi_dis(w, beta, der=0):
    r"""
    Fermi distribution.

    .. math::
       f(x) = 1/(e^x+1).

    Parameters
    ----------
    w : double
       frequency
    beta : double
       inverse temperature
    der : integer
       order of derivative

    Returns
    -------
    f : double
    """
    exponent = numpy.float128(w * beta)
    fermi = 1.0 / (numpy.exp(exponent) + 1)
    if der == 0:
        return fermi
    elif der == 1:
        return - beta * fermi ** 2 * numpy.exp(exponent)
    else:
        raise ('higher order of derivative than 1 not implemented')


def recompute_w90_input_on_different_mesh(sum_k, seedname, nk_optics, pathname='./', calc_velocity=False, calc_inverse_mass=False, oc_select='both', oc_basis='h'):
    r"""
    Recomputes dft_input objects on a finer mesh using WannierBerri and Wannier90 input.


    Parameters
    ----------
    sum_k : sum_k object
            triqs SumkDFT object
    seedname: string
              Wannier90 seedname
    nk_optics: single integer/float or three integers
               if single integer given, mesh is [nk_optics, nk_optics, nk_optics]
               elif single float given, mesh is ceiling of *sum_k.kpts * nk_optics
               elif three integers given, mesh is nk_optics
    pathname : string, optional, default='./'
               location of Wannier90 data
    calc_velocity : boolean, optional, default=False
                    whether the velocity (first derivative of H(k)) is computed
    calc_inverse_mass : boolean, optional, default=False
                        whether the inverse effective mass (second derivative of H(k)) is computed
    oc_select : string, optional, default='both'
                select contributions for optical conductivity from ['intra', 'inter', 'both']
    oc_basis : string, optional, default='h'
               gauge choice options 'h' for Hamiltonian/band and 'w' for Wannier basis

    Returns
    -------
    sum_k : sum_k object
            triqs SumkDFT object
    things_to_store : dictionary
                      dictionary of datasets to be temporarily overwritten
    """

    mpi.report('Starting Wannier interpolation...')

    BOHRTOANG = cst.physical_constants['Bohr radius'][0]/cst.angstrom
    HARTREETOEV = cst.physical_constants['Hartree energy'][0]/cst.eV
    n_inequiv_spin_blocks = sum_k.SP + 1 - sum_k.SO

    # set-up k mesh depending on input shape
    # read in transport input and some checks
    read_transport_input_from_hdf(sum_k)

    # first check for right formatting of sum_k.nk_optics
    assert len(nk_optics) in [1, 3], '"nk_optics" must be given as three integers or one float'
    if len(nk_optics) == 1:
        assert numpy.array(list(nk_optics)).dtype in (
            int, float), '"nk_optics" single value must be float or integer'
    if len(nk_optics) == 3:
        assert numpy.array(list(nk_optics)).dtype == int, '"nk_optics" mesh must be integers'
    if len(nk_optics) == 1:
        interpolate_factor = nk_optics[0]
        nk_x, nk_y, nk_z = list(
            map(lambda i: int(numpy.ceil(interpolate_factor * len(set(sum_k.kpts[:, i])))), range(3)))
    else:
        nk_x, nk_y, nk_z = nk_optics

    # check for spin calculation (not supported)
    assert sum_k.SP == 0, 'spin dependent transport calculations are not supported.'

    n_orb = numpy.max([sum_k.n_orbitals[ik][0] for ik in range(sum_k.n_k)])

    # temporarily recompute the following quantities on a different mesh
    things_to_modify = {'bz_weights': None, 'hopping': None, 'kpt_weights': None, 'kpts': None,
                        'n_k': None, 'n_orbitals': None, 'proj_mat': None, 'band_window': None, 'band_window_optics': None}
    things_to_store = dict.fromkeys(things_to_modify, None)

    # initialize variables
    n_kpts = nk_x * nk_y * nk_z
    kpts = numpy.zeros((n_kpts, 3))
    hopping = numpy.zeros((n_kpts, 1, n_orb, n_orb), dtype=complex)
    proj_mat = numpy.zeros(numpy.shape(
        hopping[:, 0, 0, 0]) + numpy.shape(sum_k.proj_mat[0, :]), dtype=complex)
    cell_volume = kpts = None
    if calc_velocity:
        velocities_k = None
    if calc_inverse_mass:
        inverse_mass = None

    if mpi.is_master_node():
        # try wannierberri import
        try:
            import wannierberri as wb
        except ImportError:
            print('ImportError: WannierBerri needs to be installed to run test "Py_w90_optics_Sr2RuO4"')
            try:
                mpi.MPI.COMM_WORLD.Abort(1)
            except:
                sys.exit()
        # initialize WannierBerri system
        shift_gamma = numpy.array([0.0, 0.0, 0.0])
        # wberri = wb.System_w90(pathname + seedname, berry=True, fft='numpy')
        # WannierBerri uses python multiprocessing which might conflict with mpi.
        # if there's a segfault, uncomment the following line
        wberri = wb.System_w90(pathname + seedname, berry=True, fft='numpy', npar=16)
        grid = wb.Grid(wberri, NKdiv=1, NKFFT=[nk_x, nk_y, nk_z])
        dataK = wb.data_K.Data_K(wberri, dK=shift_gamma, grid=grid, fftlib='numpy')

        assert dataK.HH_K.shape == hopping[:, 0, :, :].shape, 'wberri / wannier Hamiltonian has different number of orbitals than SumK object. Disentanglement is not supported as of now.'

        # read in hoppings and proj_mat
        if oc_basis == 'h':
            hopping[:, 0, range(hopping.shape[2]), range(hopping.shape[3])] = dataK.E_K
        elif oc_basis == 'w':
            hopping[:, 0, :, :] = dataK.HH_K
            fake_proj_mat = numpy.zeros(numpy.shape(dataK.UU_K), dtype=complex)
            fake_proj_mat[:, range(numpy.shape(fake_proj_mat)[1]), range(
                numpy.shape(fake_proj_mat)[2])] = 1. + 1j*0.0

        for isp in range(n_inequiv_spin_blocks):
            iorb = 0
            for icrsh in range(sum_k.n_corr_shells):
                dim = sum_k.corr_shells[icrsh]['dim']
                if oc_basis == 'h':
                    proj_mat[:, isp, icrsh, 0:dim, :] = dataK.UU_K[:, iorb:iorb+dim, :]
                elif oc_basis == 'w':
                    proj_mat[:, isp, icrsh, 0:dim, :] = fake_proj_mat[:, iorb:iorb+dim, :]
                iorb += dim

        if calc_velocity:
            # velocity: [k x n_orb x n_orb x R]
            def _commutator(A, B):
                term1 = numpy.einsum('kmo, kona -> kmna', A, B)
                term2 = numpy.einsum('kmoa, kon -> kmna', B, A)
                return term1 - term2

            # in the band basis
            # vh_alpha = Hhbar_alpha + i [Hh, Ahbar_alpha]
            if oc_basis == 'h':
                # first term
                Hhbar_alpha = dataK.Xbar('Ham', 1)
                # second term
                c_Hh_Ahbar_alpha = _commutator(hopping[:, 0, :, :], dataK.Xbar('AA'))
                velocities_k = (Hhbar_alpha + 1j * c_Hh_Ahbar_alpha) / HARTREETOEV / BOHRTOANG

                # split into diag and offdiag elements, corresponding to intra- and interband contributions
                v_diag = numpy.zeros(numpy.shape(velocities_k), dtype=complex)
                v_diag[:, range(numpy.shape(velocities_k)[1]),
                       range(numpy.shape(velocities_k)[2]), :] = velocities_k[:, range(numpy.shape(velocities_k)[1]),
                                                                              range(numpy.shape(velocities_k)[2]), :]
                v_offdiag = velocities_k.copy()
                v_offdiag[:, range(numpy.shape(velocities_k)[1]), range(
                    numpy.shape(velocities_k)[2]), :] = 0. + 1j*0.0

                if oc_select == 'intra':
                    velocities_k = v_diag
                elif oc_select == 'inter':
                    velocities_k = v_offdiag
                elif oc_select == 'both':
                    velocities_k = v_diag + v_offdiag

            # in the orbital basis
            # vw_alpha = Hw_alpha + i [Hw, Aw_alpha]
            elif oc_basis == 'w':
                # first term
                Hw_alpha_R = dataK.Ham_R.copy()
                # following three lines copied from wannierberri/data_K.py
                shape_cR = numpy.shape(dataK.cRvec_wcc)
                Hw_alpha_R = 1j * Hw_alpha_R.reshape((Hw_alpha_R.shape) + (1, )) * dataK.cRvec_wcc.reshape(
                    (shape_cR[0], shape_cR[1], dataK.system.nRvec) + (1, ) * len(Hw_alpha_R.shape[3:]) + (3, ))
                Hw_alpha = dataK.fft_R_to_k(Hw_alpha_R, hermitean=False)[dataK.select_K]
                # second term
                Aw_alpha = dataK.fft_R_to_k(dataK.AA_R, hermitean=True)
                c_Hw_Aw_alpha = _commutator(hopping[:, 0, :, :], Aw_alpha)
                velocities_k = (Hw_alpha + 1j * c_Hw_Aw_alpha) / HARTREETOEV / BOHRTOANG

        if calc_inverse_mass:
            V_dot_D = numpy.einsum('kmnab, knoab -> kmoab', dataK.Xbar('Ham', 1)
                                   [:, :, :, :, None], dataK.D_H[:, :, :, None, :])
            V_dot_D_dagger = V_dot_D.conj().transpose(0, 2, 1, 3, 4)
            V_curly = numpy.einsum('knnab -> knab', V_dot_D + V_dot_D_dagger)
            del2E_H_diag = numpy.einsum('knnab->knab', dataK.Xbar('Ham', 2)).real
            inverse_mass = del2E_H_diag + V_curly

        # read in rest from dataK
        cell_volume = dataK.cell_volume / BOHRTOANG ** 3
        kpts = dataK.kpoints_all

    # broadcast everything
    sum_k.cell_vol = mpi.bcast(cell_volume)
    kpts = mpi.bcast(kpts)
    hopping = mpi.bcast(hopping)
    proj_mat = mpi.bcast(proj_mat)
    if calc_velocity:
        velocities_k = mpi.bcast(velocities_k)
    if calc_inverse_mass:
        inverse_mass = mpi.bcast(inverse_mass)

    # write interpolated sumk quantities into "things_to_modify"
    things_to_modify['n_k'] = n_kpts
    things_to_modify['n_orbitals'] = numpy.full((n_kpts, 1), n_orb)
    for key in ['bz_weights', 'kpt_weights']:
        things_to_modify[key] = numpy.full(n_kpts, 1/n_kpts)
    n_inequiv_spin_blocks = sum_k.SP + 1 - sum_k.SO
    for key in ['band_window', 'band_window_optics']:
        things_to_modify[key] = [numpy.full((n_kpts, 2), sum_k.band_window[isp][0])
                                 for isp in range(n_inequiv_spin_blocks)]
    things_to_modify['kpts'] = kpts
    things_to_modify['hopping'] = hopping
    things_to_modify['proj_mat'] = proj_mat
    if calc_velocity:
        sum_k.velocities_k = None
        things_to_modify['velocities_k'] = velocities_k
    if calc_inverse_mass:
        sum_k.inverse_mass = None
        things_to_modify['inverse_mass'] = inverse_mass

    # now save previous sum_k instances into "things_to_store" and overwrite
    # TODO: decide whether this should be undone after the run
    for key in things_to_modify:
        things_to_store[key] = getattr(sum_k, key)
        setattr(sum_k, key, things_to_modify[key])

    # write velocities to file
    if calc_velocity:
        if mpi.is_master_node():
            ar = HDFArchive(sum_k.hdf_file, 'a')
            ar['dft_transp_input']['velocities_k'] = velocities_k
    if calc_inverse_mass:
        if mpi.is_master_node():
            ar = HDFArchive(sum_k.hdf_file, 'a')
            ar['dft_transp_input']['inverse_mass'] = inverse_mass

    return sum_k, things_to_store

# ----------------- transport -----------------------


def init_spectroscopy(sum_k, code='wien2k', w90_params={}):
    r"""
    Reads all necessary quantities for transport calculations from transport subgroup of the hdf5 archive.
    Performs checks on input. Uses interpolation if code=wannier90.

    Parameters
    ----------
    sum_k : sum_k object
            triqs SumkDFT object
    code : string
        DFT code from which velocities are being read. Options: 'wien2k', 'wannier90'
    w90_params : dictionary, optional
        additional keywords necessary in case code == 'wannier90'

    Returns
    -------
    sum_k : sum_k object
            triqs SumkDFT object, interpolated
    """

    mpi.report('Initializing optical conductivity...')
    # up and down are equivalent if SP = 0
    n_inequiv_spin_blocks = sum_k.SP + 1 - sum_k.SO

    # ----------------- set-up input from DFT -----------------------
    if code in ('wien2k', 'elk'):
        # Check if wien converter was called and read transport subgroup form
        # hdf file
        if mpi.is_master_node():
            ar = HDFArchive(sum_k.hdf_file, 'r')
            if not (sum_k.transp_data in ar):
                raise IOError(
                    "transport_distribution: No %s subgroup in hdf file found! Call convert_transp_input first." % sum_k.transp_data)
            # check if outputs file was converted
            if not ('n_symmetries' in ar['dft_misc_input']):
                raise IOError(
                    "transport_distribution: n_symmetries missing. Check if case.outputs file is present and call convert_misc_input() or convert_dft_input().")

        sum_k = read_transport_input_from_hdf(sum_k)

    elif code in ('wannier90'):
        required_entries = ['seedname', 'nk_optics']
        assert all(entry in w90_params for entry in required_entries), 'Please provide additional keywords "seedname" and "nk_optics" for "code = "wannier90""'
        # check if spin-unpolarized
        assert n_inequiv_spin_blocks == 1, "Spin-polarized optical conductivity calculations not implemented with Wannier90"

        # check some of the input
        pathname = w90_params['pathname'] if 'pathname' in w90_params else './'
        assert all(isinstance(name, str) for name in [
                   'seedname', 'pathname']), f'Check pathname {w90_params["pathname"]} and seedname {w90_params["seedname"]}'
        for file_ending in ['.wout', '_hr.dat', '.chk', '.mmn', '.eig']:
            filename = [pathname, w90_params['seedname'], file_ending]
            assert os.path.isfile(
                ''.join(filename)), f'Filename {"".join(filename)} does not exist!'
        calc_velocity = w90_params['calc_velocity'] if 'calc_velocity' in w90_params else True
        calc_inverse_mass = w90_params['calc_inverse_mass'] if 'calc_inverse_mass' in w90_params else False
        assert all(isinstance(name, bool) for name in [
                   calc_velocity, calc_inverse_mass]), f'Parameter {calc_velocity} or {calc_inverse_mass} not bool!'

        # select contributions to be used
        oc_select = w90_params['oc_select'] if 'oc_select' in w90_params else 'both'
        assert oc_select in ['intra', 'inter',
                             'both'], '"oc_select" needs to be either ["intra", "inter", "both"]'
        # gauge choice options 'h' for Hamiltonian and 'w' for Wannier
        oc_basis = w90_params['oc_basis'] if 'oc_basis' in w90_params else 'h'
        assert oc_basis in ['h', 'w'], '"oc_basis" needs to be either ["h", "w"]'
        # finally, make sure oc_select is 'both' for oc_basis = 'w'
        if oc_basis == 'w' and oc_select != 'both':
            warn(f'"oc_select" must be "both" for "oc_basis" = "w"!')
            oc_select = 'both'
        # further checks for calc_inverse_mass
        if calc_inverse_mass:
            assert oc_basis == 'h', '"calc_inverse_mass" only implemented for "oc_basis" == "h"'
            assert oc_select == 'both', '"oc_select" not implemented for "calc_inverse_mass"'
        # print some information
        mpi.report(f'{"Basis choice [h (Hamiltonian), w (Wannier)]:":<60s} {oc_basis}')
        mpi.report(f'{"Contributions from [intra(-band), inter(-band), both]:":<60s} {oc_select}')

        # recompute sum_k instances on denser grid
        sum_k, _ = recompute_w90_input_on_different_mesh(sum_k, w90_params['seedname'], nk_optics=w90_params['nk_optics'], pathname=pathname,
                                                         calc_velocity=calc_velocity, calc_inverse_mass=calc_inverse_mass, oc_select=oc_select, oc_basis=oc_basis)

    # k-dependent-projections.
    # to be checked. But this should be obsolete atm, works for both cases
    # k_dep_projection is nowhere used
    # assert sum_k.k_dep_projection == 0, "transport_distribution: k dependent projection is not implemented!"

    return sum_k

# Uses .data of only GfReFreq objects.


def transport_distribution(sum_k, beta, directions=['xx'], energy_window=None, Om_mesh=[0.0], with_Sigma=False, n_om=None, broadening=0.0, code='wien2k'):
    r"""
    Calculates the transport distribution

    .. math::
       \Gamma_{\alpha\beta}\left(\omega+\Omega/2, \omega-\Omega/2\right) = \frac{1}{V} \sum_k Tr\left(v_{k,\alpha}A_{k}(\omega+\Omega/2)v_{k,\beta}A_{k}\left(\omega-\Omega/2\right)\right)

    in the direction :math:`\alpha\beta`. The velocities :math:`v_{k}` are read from the transport subgroup of the hdf5 archive.

    Parameters
    ----------
    sum_k : sum_k object
            triqs SumkDFT object
    beta : double
        Inverse temperature :math:`\beta`.
    directions : list of string, optional
        :math:`\alpha\beta` e.g.: ['xx','yy','zz','xy','xz','yz'].
    energy_window : list of double, optional
        Specifies the upper and lower limit of the frequency integration for :math:`\Omega=0.0`. The window is automatically enlarged by the largest :math:`\Omega` value,
        hence the integration is performed in the interval [energy_window[0]-max(Om_mesh), energy_window[1]+max(Om_mesh)].
    Om_mesh : list of double, optional
        :math:`\Omega` frequency mesh of the optical conductivity. For the conductivity and the Seebeck coefficient :math:`\Omega=0.0` has to be
        part of the mesh. In the current version Om_mesh is repined to the mesh provided by the self-energy! The actual mesh is printed on the screen and given as output.
    with_Sigma : boolean, optional
        Determines whether the calculation is performed with or without self energy. If this parameter is set to False the self energy is set to zero (i.e. the DFT band
        structure :math:`A(k,\omega)` is used). Note: For with_Sigma=False it is necessary to specify the parameters energy_window, n_om and broadening.
    n_om : integer, optional
        Number of equidistant frequency points in the interval [energy_window[0]-max(Om_mesh), energy_window[1]+max(Om_mesh)]. This parameters is only used if
        with_Sigma = False.
    broadening : double, optional
        Lorentzian broadening. It is necessary to specify the boradening if with_Sigma = False, otherwise this parameter can be set to 0.0.
    code : string
        DFT code from which velocities are being read. Options: 'wien2k', 'wannier90'

    Returns
    -------
    Gamma_w : dictionary of double matrices
              transport distribution function in each direction, frequency given by Om_mesh_out and omega
    omega : list of double
            omega vector
    Om_mesh_out : list of double
                  frequency mesh of the optical conductivity recomputed on the mesh provided by the self energy
    """

    mpi.report('Computing transport distribution...')

    n_inequiv_spin_blocks = sum_k.SP + 1 - sum_k.SO
    # up and down are equivalent if SP = 0

    # positive om_mesh
    assert all(
        Om >= 0.0 for Om in Om_mesh), "transport_distribution: Om_mesh should not contain negative values!"
    # Check if energy_window is sufficiently large and correct
    if (energy_window[0] >= energy_window[1] or energy_window[0] >= 0 or energy_window[1] <= 0):
        assert 0, "transport_distribution: energy_window wrong!"

    if (abs(fermi_dis(energy_window[0], beta) * fermi_dis(-energy_window[0], beta)) > 1e-5
            or abs(fermi_dis(energy_window[1], beta) * fermi_dis(-energy_window[1], beta)) > 1e-5):
        mpi.report(
            "\n####################################################################")
        mpi.report(
            "transport_distribution: WARNING - energy window might be too narrow!")
        mpi.report(
            "####################################################################\n")

    # ----------------- calculate A(k,w) -----------------------

    # Define mesh for Green's function and in the specified energy window
    if (with_Sigma == True):
        omega = numpy.array([round(x.real, 12)
                             for x in sum_k.Sigma_imp[0].mesh])
        mesh = None
        mu = sum_k.chemical_potential
        n_om = len(omega)
        mpi.report("Using omega mesh provided by Sigma!")

        if energy_window:
            # Find according window in Sigma mesh
            ioffset = numpy.sum(
                omega < energy_window[0] - max(Om_mesh))
            omega = omega[numpy.logical_and(
                omega >= energy_window[0] - max(Om_mesh), omega <= energy_window[1] + max(Om_mesh))]
            n_om = len(omega)

            # Truncate Sigma to given omega window
            # In the future there should be an option in gf to manipulate the mesh (e.g. truncate) directly.
            # For now we stick with this:
            for icrsh in range(sum_k.n_corr_shells):
                Sigma_save = sum_k.Sigma_imp[icrsh].copy()
                spn = sum_k.spin_block_names[sum_k.corr_shells[icrsh]['SO']]
                def glist(): return [GfReFreq(target_shape=(block_dim, block_dim), window=(omega[
                    0], omega[-1]), n_points=n_om) for block, block_dim in sum_k.gf_struct_sumk[icrsh]]
                sum_k.Sigma_imp[icrsh] = BlockGf(
                    name_list=spn, block_list=glist(), make_copies=False)
                for i, g in sum_k.Sigma_imp[icrsh]:
                    for iL in g.indices[0]:
                        for iR in g.indices[0]:
                            for iom in range(n_om):
                                g.data[iom, int(iL), int(iR)] = Sigma_save[
                                    i].data[ioffset + iom, int(iL), int(iR)]
    else:
        assert n_om is not None, "transport_distribution: Number of omega points (n_om) needed to calculate transport distribution!"
        assert energy_window is not None, "transport_distribution: Energy window needed to calculate transport distribution!"
        assert broadening != 0.0 and broadening is not None, "transport_distribution: Broadening necessary to calculate transport distribution!"
        omega = numpy.linspace(
            energy_window[0] - max(Om_mesh), energy_window[1] + max(Om_mesh), n_om)
        mesh = MeshReFreq(energy_window[0] -
                          max(Om_mesh), energy_window[1] + max(Om_mesh), n_om)
        mu = 0.0

    dir_to_int = {'x': 0, 'y': 1, 'z': 2}

    # Define mesh for optic conductivity
    d_omega = round(numpy.abs(omega[0] - omega[1]), 12)
    iOm_mesh = numpy.array([round((Om / d_omega), 0) for Om in Om_mesh])
    temp_Om_mesh = iOm_mesh * d_omega

    if mpi.is_master_node():
        print("Chemical potential: ", mu)
        print("Using n_om = %s points in the energy_window [%s,%s]" % (
            n_om, omega[0], omega[-1]), end=' ')
        print("where the omega vector is:")
        print(omega)
        print("Calculation requested for Omega mesh:   ", numpy.array(Om_mesh))
        print("Omega mesh automatically repined to:  ", temp_Om_mesh)

    Gamma_w = {direction: numpy.zeros((len(temp_Om_mesh), n_om),
                                      dtype=numpy.float_) for direction in directions}

    # Sum over all k-points
    ikarray = numpy.array(list(range(sum_k.n_k)))
    for ik in mpi.slice_array(ikarray):
        # Calculate G_w  for ik and initialize A_kw
        G_w = sum_k.lattice_gf(ik, mu, broadening=broadening, mesh=mesh, with_Sigma=with_Sigma)
        A_kw = [numpy.zeros((sum_k.n_orbitals[ik][isp], sum_k.n_orbitals[ik][isp], n_om), dtype=numpy.complex_)
                for isp in range(n_inequiv_spin_blocks)]

        for isp in range(n_inequiv_spin_blocks):
            # copy data from G_w (swapaxes is used to have omega in the 3rd
            # dimension)
            A_kw[isp] = copy.deepcopy(G_w[sum_k.spin_block_names[sum_k.SO][
                                      isp]].data.swapaxes(0, 1).swapaxes(1, 2))
            # calculate A(k,w) for each frequency
            for iw in range(n_om):
                A_kw[isp][:, :, iw] = -1.0 / (2.0 * numpy.pi * 1j) * (
                    A_kw[isp][:, :, iw] - numpy.conjugate(numpy.transpose(A_kw[isp][:, :, iw])))
            # Akw_write[ik] = A_kw[isp].copy() * sum_k.bz_weights[ik]

            b_min = max(sum_k.band_window[isp][
                        ik, 0], sum_k.band_window_optics[isp][ik, 0])
            b_max = min(sum_k.band_window[isp][
                        ik, 1], sum_k.band_window_optics[isp][ik, 1])
            A_i = slice(
                b_min - sum_k.band_window[isp][ik, 0], b_max - sum_k.band_window[isp][ik, 0] + 1)
            v_i = slice(b_min - sum_k.band_window_optics[isp][
                        ik, 0], b_max - sum_k.band_window_optics[isp][ik, 0] + 1)

            # loop over all symmetries
            for R in sum_k.rot_symmetries:
                # get transformed velocity under symmetry R
                if code in ('wien2k'):
                    vel_R = copy.deepcopy(sum_k.velocities_k[isp][ik])
                elif code in ('wannier90'):
                    vel_R = copy.deepcopy(sum_k.velocities_k[ik])
                for nu1 in range(sum_k.band_window_optics[isp][ik, 1] - sum_k.band_window_optics[isp][ik, 0] + 1):
                    for nu2 in range(sum_k.band_window_optics[isp][ik, 1] - sum_k.band_window_optics[isp][ik, 0] + 1):
                        vel_R[nu1][nu2][:] = numpy.dot(
                            R, vel_R[nu1][nu2][:])

                # calculate Gamma_w for each direction from the velocities
                # vel_R and the spectral function A_kw
                for direction in directions:
                    for iw in range(n_om):
                        for iq in range(len(temp_Om_mesh)):
                            if (iw + iOm_mesh[iq] >= n_om or omega[iw] < -temp_Om_mesh[iq] + energy_window[0] or omega[iw] > temp_Om_mesh[iq] + energy_window[1]):
                                continue

                            Gamma_w[direction][iq, iw] += (numpy.dot(numpy.dot(numpy.dot(vel_R[v_i, v_i, dir_to_int[direction[0]]], A_kw[isp][A_i, A_i, int(iw + iOm_mesh[iq])]),
                                                                               vel_R[v_i, v_i, dir_to_int[direction[1]]]), A_kw[isp][A_i, A_i, iw]).trace().real * sum_k.bz_weights[ik])

    for direction in directions:
        Gamma_w[direction] = (mpi.all_reduce(mpi.world, Gamma_w[direction],
                              lambda x, y: x + y) / sum_k.cell_vol / sum_k.n_symmetries)

    return Gamma_w, omega, temp_Om_mesh


def transport_function(beta, directions, hopping, velocities, energy_window, n_om, rot_symmetries):
    r"""
    Calculates the transport function

    .. math::
       \Phi_\alpha\beta(\omega) = \sum_k v_{k,\alpha} v_{k,\beta} \delta(\omega-\varepsilon)

    in the direction :math:`\alpha\beta`.

    Parameters
    ----------
    beta : double
        Inverse temperature :math:`\beta`.
    directions : list of string, optional
        :math:`\alpha\beta` e.g.: ['xx','yy','zz','xy','xz','yz'].
    hopping : double array
        Hamiltonian in band basis :math:`\epsilon(k)`
    veolcities : complex array
        matrix elements derivative of Hamiltonian :math:`\frac{d\epsilon(k)}{dk}`
    energy_window : list of double
        Specifies the upper and lower limit of the frequency integration for :math:`\Omega=0.0`. The window is automatically enlarged by the largest :math:`\Omega` value,
        hence the integration is performed in the interval [energy_window[0]-max(Om_mesh), energy_window[1]+max(Om_mesh)].
    n_om : integer
        Number of equidistant frequency points in the interval [energy_window[0]-max(Om_mesh), energy_window[1]+max(Om_mesh)]. This parameters is only used if
        with_Sigma = False.
    rot_symmetries :  list of 3 x 3 matrices
        rotational symmetries to restore the full FBZ

    Returns
    -------
    transp_func : dictionary of double array
              transport function in each direction, frequencies given by energy_window
    """

    mpi.report('Computing transport function...')

    # check that velocities are computed on the FBZ
    assert numpy.shape(rot_symmetries)[
        0] == 1, 'Using symmetries currently not implemented for transport function.'

    dir_to_int = {'x': 0, 'y': 1, 'z': 2}

    tol = 1/beta
    orb_1, orb_2 = velocities.shape[1:3]
    ws = numpy.linspace(energy_window[0], energy_window[1], n_om)
    transp_func = {direction: numpy.zeros(shape=(ws.shape[0])) for direction in directions}

    for ct, w in enumerate(ws):
        idx = numpy.where(numpy.abs(hopping[:, 0, range(orb_1), range(orb_2)].real - w) <= tol)
        fermi_wg = fermi_dis(hopping[:, 0, range(orb_1), range(orb_2)]
                             [idx].real - w, beta, 1)/fermi_dis(0., beta, 1)
        for direction in directions:
            dir_a, dir_b = [dir_to_int[x] for x in direction]
            matrix_product = numpy.einsum(
                'kmn, kno -> kmo', velocities[:, :, :, dir_a], velocities[:, :, :, dir_b])
            transp_func[direction][ct] = numpy.sum(
                fermi_wg * matrix_product[:, range(orb_1), range(orb_2)][idx]).real

    return transp_func


def transport_coefficient(Gamma_w, omega, Om_mesh, spin_polarization, direction, iq, n, beta, method=None):
    r"""
    Calculates the transport coefficient A_n in a given direction for a given :math:`\Omega`. The required members (Gamma_w, directions, Om_mesh) have to be obtained first
    by calling the function :meth:`transport_distribution <dft.sumk_dft_tools.SumkDFTTools.transport_distribution>`. For n>0 A is set to NaN if :math:`\Omega` is not 0.0.

    Parameters
    ----------
    Gamma_w : dictionary of double matrices
              transport distribution function in each direction, frequency given by Om_mesh_out and omega
    omega : list of double
            omega vector
    Om_mesh : list of double
              frequency mesh of the optical conductivity recomputed on the mesh provided by the self energy
    spin_polarization : integer
                        Boolean-type integer whether system is spin-polarized or not
    direction : string
       :math:`\alpha\beta` e.g.: 'xx','yy','zz','xy','xz','yz'.
    iq : integer
        Index of :math:`\Omega` point in the member Om_mesh.
    n : integer
        Number of the desired moment of the transport distribution.
    beta : double
        Inverse temperature :math:`\beta`.
    method : string
        Integration method: cubic spline and scipy.integrate.quad ('quad'), simpson rule ('simps'), trapezoidal rule ('trapz'), rectangular integration (otherwise)
        Note that the sampling points of the the self-energy are used!

    Returns
    -------
    A : double
        Transport coefficient.
    """

    from scipy.interpolate import interp1d
    from scipy.integrate import simpson, quad

    if not (mpi.is_master_node()):
        return

    if (Om_mesh[iq] == 0.0 or n == 0.0):
        A = 0.0
        # setup the integrand
        if (Om_mesh[iq] == 0.0):
            A_int = Gamma_w[direction][iq] * \
                (fermi_dis(omega, beta) * fermi_dis(-omega, beta)) * (omega * beta)**n
        elif (n == 0.0):
            A_int = Gamma_w[direction][iq] * (fermi_dis(omega, beta) -
                                              fermi_dis(omega + Om_mesh[iq], beta)) / (Om_mesh[iq] * beta)

        # w-integration
        if method == 'quad':
            # quad on interpolated w-points with cubic spline
            A_int_interp = interp1d(omega, A_int, kind='cubic')
            A = quad(A_int_interp, min(omega), max(omega),
                     epsabs=1.0e-12, epsrel=1.0e-12, limit=500)
            A = A[0]
        elif method == 'simps':
            # simpson rule for w-grid
            A = simpson(A_int, omega)
        elif method == 'trapz':
            # trapezoidal rule for w-grid
            A = numpy.trapz(A_int, omega)
        else:
            # rectangular integration for w-grid (orignal implementation)
            d_w = omega[1] - omega[0]
            for iw in range(Gamma_w[direction].shape[1]):
                A += A_int[iw] * d_w
        A = A * numpy.pi * (2.0 - spin_polarization)
    else:
        A = numpy.nan
    return A


def conductivity_and_seebeck(Gamma_w, omega, Om_mesh, SP, directions, beta, method=None):
    r"""
    Calculates the Seebeck coefficient and the optical conductivity by calling
    :meth:`transport_coefficient <dft.sumk_dft_tools.SumkDFTTools.transport_coefficient>`.
    The required members (Gamma_w, directions, Om_mesh) have to be obtained first by calling the function
    :meth:`transport_distribution <dft.sumk_dft_tools.SumkDFTTools.transport_distribution>`.

    Parameters
    ----------
    Gamma_w : dictionary of double matrices
              transport distribution function in each direction, frequency given by Om_mesh_out and omega
    omega : list of double
            omega vector
    Om_mesh : list of double
              frequency mesh of the optical conductivity recomputed on the mesh provided by the self energy
    spin_polarization : integer
                        Boolean-type integer whether system is spin-polarized or not
    directions : list of string, optional
        :math:`\alpha\beta` e.g.: ['xx','yy','zz','xy','xz','yz'].
    beta : double
        Inverse temperature :math:`\beta`.
    method : string
        Integration method: cubic spline and scipy.integrate.quad ('quad'), simpson rule ('simps'), trapezoidal rule ('trapz'), rectangular integration (otherwise)
        Note that the sampling points of the the self-energy are used!

    Returns
    -------
    optic_cond : dictionary of double vectors
        Optical conductivity in each direction and frequency given by Om_mesh.

    seebeck : dictionary of double
        Seebeck coefficient in each direction. If zero is not present in Om_mesh the Seebeck coefficient is set to NaN.

    kappa : dictionary of double.
        thermal conductivity in each direction. If zero is not present in Om_mesh the thermal conductivity is set to NaN
    """

    mpi.report('Computing optical conductivity and kinetic coefficients...')

    if not (mpi.is_master_node()):
        return

    n_q = Gamma_w[directions[0]].shape[0]

    # initialization
    A0 = {direction: numpy.full((n_q,), numpy.nan) for direction in directions}
    A1 = {direction: numpy.full((n_q,), numpy.nan) for direction in directions}
    A2 = {direction: numpy.full((n_q,), numpy.nan) for direction in directions}
    optic_cond = {direction: numpy.full((n_q,), numpy.nan) for direction in directions}
    seebeck = {direction: numpy.nan for direction in directions}
    kappa = {direction: numpy.nan for direction in directions}

    for direction in directions:
        for iq in range(n_q):
            A0[direction][iq] = transport_coefficient(
                Gamma_w, omega, Om_mesh, SP, direction, iq=iq, n=0, beta=beta, method=method)
            A1[direction][iq] = transport_coefficient(
                Gamma_w, omega, Om_mesh, SP, direction, iq=iq, n=1, beta=beta, method=method)
            A2[direction][iq] = transport_coefficient(
                Gamma_w, omega, Om_mesh, SP, direction, iq=iq, n=2, beta=beta, method=method)
            print("A_0 in direction %s for Omega = %.2f    %e a.u." %
                  (direction, Om_mesh[iq], A0[direction][iq]))
            print("A_1 in direction %s for Omega = %.2f    %e a.u." %
                  (direction, Om_mesh[iq], A1[direction][iq]))
            print("A_2 in direction %s for Omega = %.2f    %e a.u." %
                  (direction, Om_mesh[iq], A2[direction][iq]))
            if ~numpy.isnan(A1[direction][iq]):
                # Seebeck and kappa are overwritten if there is more than one Omega =
                # 0 in Om_mesh
                seebeck[direction] = - A1[direction][iq] / A0[direction][iq] * 86.17
                kappa[direction] = A2[direction][iq] - \
                    A1[direction][iq]*A1[direction][iq]/A0[direction][iq]
                kappa[direction] *= 293178.0

        # factor for optical conductivity: hbar * velocity_Hartree_to_SI * volume_Hartree_to_SI * m_to_cm * 10^-4 final unit
        convert_to_SI = cst.hbar * (cst.c * cst.fine_structure) ** 2 * \
            (1/cst.physical_constants['Bohr radius'][0]) ** 3 * 1e-6
        optic_cond[direction] = beta * convert_to_SI * A0[direction]
        for iq in range(n_q):
            print("Conductivity in direction %s for Omega = %.2f       %f  x 10^4 Ohm^-1 cm^-1" %
                  (direction, Om_mesh[iq], optic_cond[direction][iq]))
            if not (numpy.isnan(A1[direction][iq])):
                print("Seebeck in direction      %s for Omega = 0.00      %f  x 10^(-6) V/K" %
                      (direction, seebeck[direction]))
                print("kappa in direction      %s for Omega = 0.00      %f  W/(m * K)" %
                      (direction, kappa[direction]))

    return optic_cond, seebeck, kappa
