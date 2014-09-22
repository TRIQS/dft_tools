from math import sqrt
from scipy.misc import factorial as fact
from itertools import product
import numpy as np


# The interaction matrix in desired basis
# U^{spherical}_{m1 m2 m3 m4} = \sum_{k=0}^{2l} F_k angular_matrix_element(l, k, m1, m2, m3, m4)
def U_matrix(l, radial_integrals=None, U_int=None, J_hund=None, basis="spherical", T=None):
    """Calculate U matrix being given either radial_integrals or U_int and J_hund.
       l = angular momentum of shell being treated (l=2 for d shell, l=3 for f shell)
       radial_integrals = [F0,F2,F4,..] (default None)
       U_int, J_hund = values to use to compute radial_integrals (default None),
       basis = "spherical", "cubic", or "other",
       T = transformation matrix from spherical to desired basis, if basis='other' (default None)"""

    # Check all necessary information is present and consistent
    if radial_integrals == None and (U_int == None and J_hund == None):
        raise ValueError("U_matrix: provide either the radial_integrals or U_int and J_hund.")
    if radial_integrals == None and (U_int != None and J_hund != None):
        radial_integrals = U_J_to_radial_integrals(l, U_int, J_hund)
    if radial_integrals != None and (U_int != None and J_hund != None):
        if len(radial_integrals)-1 != l:
            raise ValueError("U_matrix: inconsistency in l and number of radial_integrals provided.")
        if (radial_integrals - U_J_to_radial_integrals(l, U_int, J_hund)).any() != 0.0:
            print "Warning: U_matrix: radial_integrals provided do not match U_int and J_hund. Using radial_integrals to calculate spherical U_matrix."

    # Full interaction matrix
    # Basis of spherical harmonics Y_{-2}, Y_{-1}, Y_{0}, Y_{1}, Y_{2}
    # U^{spherical}_{m1 m2 m3 m4} = \sum_{k=0}^{2l} F_k angular_matrix_element(l, k, m1, m2, m3, m4)
    U_matrix = np.zeros((2*l+1,2*l+1,2*l+1,2*l+1),dtype=float)

    m_range = range(-l,l+1)
    for n, F in enumerate(radial_integrals):
        k = 2*n
        for m1, m2, m3, m4 in product(m_range,m_range,m_range,m_range):
            U_matrix[m1+l,m2+l,m3+l,m4+l] += F * angular_matrix_element(l,k,m1,m2,m3,m4)

    # Transform from spherical basis if needed
    if basis == "cubic": T = spherical_to_cubic(l)
    if basis == "other" and T == None:
        raise ValueError("U_matrix: provide T for other bases.")
    if T != None: transform_U_matrix(U_matrix, T)

    return U_matrix

# Convert full 4-index U matrix to 2-index density-density form
def reduce_4index_to_2index(U_4index):
    """Reduces the four-index matrix to two-index matrices."""

    size = len(U_4index) # 2l+1
    U  = np.zeros((size,size),dtype=float)      # matrix for same spin
    Uprime = np.zeros((size,size),dtype=float)  # matrix for opposite spin

    m_range = range(size)
    for m,mp in product(m_range,m_range):
        U[m,mp]  = U_4index[m,mp,m,mp].real - U_4index[m,mp,mp,m].real
        Uprime[m,mp] = U_4index[m,mp,m,mp].real

    return U, Uprime

# Construct the 2-index matrices for the density-density form
def U_matrix_kanamori(n_orb, U_int=None, J_hund=None):
    """Calculate the Kanamori U and Uprime matrices."""

    U  = np.zeros((n_orb,n_orb),dtype=float)      # matrix for same spin
    Uprime = np.zeros((n_orb,n_orb),dtype=float)  # matrix for opposite spin

    m_range = range(n_orb)
    for m,mp in product(m_range,m_range):
        if m == mp:
            Uprime[m,mp] = U_int
        else:
            U[m,mp]  = U_int - 3.0*J_hund
            Uprime[m,mp] = U_int - 2.0*J_hund

    return U, Uprime

# Get t2g or eg components
def t2g_submatrix(U):
    """Return only the t2g part of the full d-manifold two- or four-index U matrix."""
    return subarray(U, len(U.shape)*[(0,1,3)])

def eg_submatrix(U):
    """Return only the eg part of the full d-manifold two- or four-index U matrix."""
    return subarray(U, len(U.shape)*[(2,4)])

# Transform the interaction matrix into another basis
def transform_U_matrix(U_matrix, T):
    """Transform the interaction matrix into another basis by applying matrix T."""
    return np.einsum("ij,kl,jlmo,mn,op",np.conj(T),np.conj(T),U_matrix,np.transpose(T),np.transpose(T))

# Rotation matrices: complex harmonics to cubic harmonics
# Complex harmonics basis: ..., Y_{-2}, Y_{-1}, Y_{0}, Y_{1}, Y_{2}, ...
def spherical_to_cubic(l):
    """Returns the spherical harmonics to cubic harmonics rotation matrix."""
    size = 2*l+1
    T = np.zeros((size,size),dtype=complex)
    if l == 0:
        cubic_names = ("s")
    elif l == 1:
        cubic_names = ("x","y","z")
        T[0,0] = 1.0/sqrt(2);   T[0,2] = -1.0/sqrt(2)
        T[1,0] = 1j/sqrt(2);    T[1,2] = 1j/sqrt(2)
        T[2,1] = 1.0
    elif l == 2:
        cubic_names = ("xy","yz","z^2","xz","x^2-y^2")
        T[0,0] = 1j/sqrt(2);    T[0,4] = -1j/sqrt(2)
        T[1,1] = 1j/sqrt(2);    T[1,3] = 1j/sqrt(2)
        T[2,2] = 1.0
        T[3,1] = 1.0/sqrt(2);   T[3,3] = -1.0/sqrt(2)
        T[4,0] = 1.0/sqrt(2);   T[4,4] = 1.0/sqrt(2)
    elif l == 3:
        cubic_names = ("x(x^2-3y^2)","z(x^2-y^2)","xz^2","z^3","yz^2","xyz","y(3x^2-y^2)")
        T[0,0] = 1.0/sqrt(2);    T[0,6] = -1.0/sqrt(2)
        T[1,1] = 1.0/sqrt(2);    T[1,5] = 1.0/sqrt(2)
        T[2,2] = 1.0/sqrt(2);    T[2,4] = -1.0/sqrt(2)
        T[3,5] = 1.0
        T[4,2] = 1j/sqrt(2);   T[4,4] = 1j/sqrt(2)
        T[5,1] = 1j/sqrt(2);   T[5,5] = -1j/sqrt(2)
        T[6,0] = 1j/sqrt(2);   T[6,6] = 1j/sqrt(2)
    else: raise ValueError("spherical_to_cubic: implemented only for l=0,1,2,3")

    return T

# Names of cubic harmonics
def cubic_names(l):
    if l == 0 or l == 's':
        return ("s")
    elif l == 1 or l == 'p':
        return ("x","y","z")
    elif l == 2 or l == 'd':
        return ("xy","yz","z^2","xz","x^2-y^2")
    elif l == 't2g':
        return ("xy","yz","xz")
    elif l == 'eg':
        return ("z^2","x^2-y^2")
    elif l == 3 or l == 'f':
        return ("x(x^2-3y^2)","z(x^2-y^2)","xz^2","z^3","yz^2","xyz","y(3x^2-y^2)")
    else: raise ValueError("cubic_names: implemented only for l=0,1,2,3")

# Convert U,J -> radial integrals F_k
def U_J_to_radial_integrals(l, U_int, J_hund):
    """Determines the radial integrals F_k from U_int and J_hund."""

    F = np.zeros((l+1),dtype=float)
    if l == 2:
        F[0] = U_int
        F[1] = J_hund * 14.0 / (1.0 + 0.63)
        F[2] = 0.630 * F[1]
    elif l == 3:
        F[0] = U_int
        F[1] = 6435.0 * J_hund / (286.0 + 195.0 * 451.0 / 675.0 + 250.0 * 1001.0 / 2025.0)
        F[2] = 451.0 * F[1] / 675.0
        F[3] = 1001.0 * F[1] / 2025.0
    else: raise ValueError("U_J_to_radial_integrals: implemented only for l=2,3")

    return F

# Angular matrix elements of particle-particle interaction
# (2l+1)^2 ((l 0) (k 0) (l 0))^2 \sum_{q=-k}^{k} (-1)^{m1+m2+q} ((l -m1) (k q) (l m3)) ((l -m2) (k -q) (l m4))
def angular_matrix_element(l, k, m1, m2, m3, m4):
    result = 0
    for q in range(-k,k+1):
        result += three_j_symbol((l,-m1),(k,q),(l,m3))*three_j_symbol((l,-m2),(k,-q),(l,m4))*(-1.0 if (m1+q+m2) % 2 else 1.0)
    result *= (2*l+1)**2 * (three_j_symbol((l,0),(k,0),(l,0))**2)
    return result

# Wigner 3-j symbols
# ((j1 m1) (j2 m2) (j3 m3))
def three_j_symbol(jm1, jm2, jm3):
    j1, m1 = jm1
    j2, m2 = jm2
    j3, m3 = jm3

    if (m1+m2+m3 != 0 or
        m1 < -j1 or m1 > j1 or
        m2 < -j2 or m2 > j2 or
        m3 < -j3 or m3 > j3 or
        j3 > j1 + j2 or
        j3 < abs(j1-j2)):
        return .0

    result = -1.0 if (j1-j2-m3) % 2 else 1.0
    result *= sqrt(fact(j1+j2-j3)*fact(j1-j2+j3)*fact(-j1+j2+j3)/fact(j1+j2+j3+1))
    result *= sqrt(fact(j1-m1)*fact(j1+m1)*fact(j2-m2)*fact(j2+m2)*fact(j3-m3)*fact(j3+m3))

    t_min = max(j2-j3-m1,j1-j3+m2,0)
    t_max = min(j1-m1,j2+m2,j1+j2-j3)

    t_sum = 0
    for t in range(t_min,t_max+1):
        t_sum += (-1.0 if t % 2 else 1.0)/(fact(t)*fact(j3-j2+m1+t)*fact(j3-j1-m2+t)*fact(j1+j2-j3-t)*fact(j1-m1-t)*fact(j2+m2-t))

    result *= t_sum
    return result

# Clebsch-Gordan coefficients
# < j1 m1 j2 m2 | j3 m3 > = (-1)^{j1-j2+m3} \sqrt{2j3+1} ((j1 m1) (j2 m2) (j3 -m3))
def clebsch_gordan(jm1, jm2, jm3):
    norm = sqrt(2*jm3[0]+1)*(-1 if jm1[0]-jm2[0]+jm3[1] % 2 else 1)
    return norm*three_j_symbol(jm1,jm2,(jm3[0],-jm3[1]))

# Create subarray containing columns in idxlist
# e.g. idxlist = [(0),(2,3),(0,1,2,3)] gives
#  column 0 for 1st dim,
#  columns 2 and 3 for 2nd dim,
#  columns 0,1,2 and 3 for 3rd dim.
#def subarray(a,idxlist,n=len(a.shape)-1) :
def subarray(a,idxlist,n=None) :
    if n == None: n = len(a.shape)-1
    sa = a[tuple(slice(x) for x in a.shape[:n]) + (idxlist[n],)]
    return subarray(sa,idxlist, n-1) if n > 0 else sa
