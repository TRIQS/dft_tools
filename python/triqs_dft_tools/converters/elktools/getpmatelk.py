
##########################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2025 by N. Wentzell
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

import numpy as np
import os

def getpmatelk(ik, nstsv, vkl, filename='PMAT.OUT'):
    """
    Pure Python implementation to read momentum matrix elements from Elk's PMAT.OUT file.

    This replaces the previous f2py compiled Fortran module.

    Parameters
    ----------
    ik : int
        K-point index (1-based, as in Fortran)
    nstsv : int
        Number of states (second variation)
    vkl : array_like
        K-point lattice coordinates (3 floats)
    filename : str, optional
        Path to PMAT.OUT file (default: 'PMAT.OUT')

    Returns
    -------
    pmat : ndarray
        Momentum matrix elements with shape (nstsv, nstsv, 3), dtype=complex128
    """

    if not os.path.exists(filename):
        raise IOError(f"Error(getpmat): unable to read from {filename}")

    # Tolerance for k-point comparison (from Elk's epslat)
    epslat = 1.0e-6

    # Data types matching Fortran:
    # real(8) -> float64 (8 bytes)
    # integer -> int32 (4 bytes, standard Fortran default integer)
    #            Note: If Elk is compiled with -fdefault-integer-8, use int64
    # complex(8) -> complex128 (16 bytes: 8-byte real + 8-byte imag)

    # Calculate record length in bytes
    # Record contains: vkl(3) + nstsv + pmat(nstsv, nstsv, 3)
    vkl_bytes = 3 * 8  # 3 real(8)
    nstsv_bytes = 4    # 1 integer
    pmat_bytes = nstsv * nstsv * 3 * 16  # nstsv*nstsv*3 complex(8)
    recl = vkl_bytes + nstsv_bytes + pmat_bytes

    # Read the record at position ik (1-based indexing)
    # In direct access, record ik starts at byte position (ik-1) * recl
    try:
        with open(filename, 'rb') as f:
            # Seek to the correct record
            f.seek((ik - 1) * recl)

            # Read vkl (3 float64)
            vkl_read = np.fromfile(f, dtype=np.float64, count=3)
            if vkl_read.size != 3:
                raise IOError(f"Error(getpmat): unable to read from {filename}")

            # Read nstsv (1 int32)
            nstsv_array = np.fromfile(f, dtype=np.int32, count=1)
            if nstsv_array.size != 1:
                raise IOError(f"Error(getpmat): unable to read from {filename}")
            nstsv_read = nstsv_array[0]

            # Read pmat (nstsv*nstsv*3 complex128)
            # Fortran stores arrays in column-major order
            # pmat(nstsv, nstsv, 3) in Fortran means first index varies fastest
            pmat_flat = np.fromfile(f, dtype=np.complex128, count=nstsv*nstsv*3)
            if pmat_flat.size != nstsv*nstsv*3:
                raise IOError(f"Error(getpmat): unable to read from {filename}")

            # Reshape directly to (nstsv, nstsv, 3) with Fortran order
            pmat = pmat_flat.reshape((nstsv, nstsv, 3), order='F')
    except (OSError, IOError) as e:
        raise IOError(f"Error(getpmat): unable to read from {filename}") from e

    # Validate k-point
    vkl = np.asarray(vkl)
    t1 = np.sum(np.abs(vkl - vkl_read))
    if t1 > epslat:
        raise ValueError(
            f"Error(getpmat): differing vectors for k-point {ik}\n"
            f"  current  : {vkl}\n"
            f"  PMAT.OUT : {vkl_read}"
        )

    # Validate nstsv
    if nstsv != nstsv_read:
        raise ValueError(
            f"Error(getpmat): differing nstsv for k-point {ik}\n"
            f"  current  : {nstsv}\n"
            f"  PMAT.OUT : {nstsv_read}"
        )

    return pmat
