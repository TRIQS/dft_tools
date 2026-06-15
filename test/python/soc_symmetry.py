################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
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
################################################################################

# Regression test for the spin-orbit symmetrization (issue #148).
#
# soc_symmetry.ref.h5 holds the correlated-shell symmetry data (dmftproj
# case.symqmc) of a full-BZ spin-orbit + spin-polarized Wien2k calculation for a
# cubic fluorite-type cell with two symmetry-equivalent correlated atoms. Its
# magnetic point group has 16 operations, 8 of them time-reversal (the path whose
# phase handling #148 questions).
#
# The symmetry operations must form a (anti)unitary representation, so the
# group-average symmetrizer is a projector: applied to an already-symmetric
# matrix it must not change its eigenvalues. A sign or structure error in the
# spinor matrices breaks this and shifts the eigenvalues, which is the failure
# reported in #148.

import numpy as np
from triqs_dft_tools.symmetry import Symmetry

symm = Symmetry('soc_symmetry.ref.h5', subgroup='dft_symmcorr_input')
assert sum(int(t) for t in symm.time_inv) > 0, "fixture should contain time_inv operations"

dim = symm.mat[0][0].shape[0]
rng = np.random.RandomState(2148)
random_herm = []
for _ in range(symm.n_orbits):
    a = rng.randn(dim, dim) + 1j * rng.randn(dim, dim)
    random_herm.append(0.5 * (a + a.conj().transpose()))

# one symmetrization lands in the symmetric subspace; a second must be a no-op.
# case.symqmc stores the matrices as text with ~6 decimals, so the projector
# identity holds only to that precision; a sign/structure error breaks it by O(1).
tol = 1e-5
symmetric = symm.symmetrize(random_herm)
twice = symm.symmetrize(symmetric)

for orb in range(symm.n_orbits):
    assert np.allclose(symmetric[orb], twice[orb], atol=tol), \
        "symmetrizer is not idempotent (orbit %d): the spinor symmetry matrices " \
        "do not form a valid representation" % orb
    e1 = np.sort(np.linalg.eigvalsh(symmetric[orb]))
    e2 = np.sort(np.linalg.eigvalsh(twice[orb]))
    assert np.allclose(e1, e2, atol=tol), \
        "symmetrization changed the eigenvalues of an already-symmetric matrix " \
        "(orbit %d)" % orb
