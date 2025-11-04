import numpy as np

_TOL = 1e-8

# Ensure strictly positive imaginary part with minimal scale
_regularize = lambda z : 1j * max(float(z), 1.0e-20) / 100.0

def _F(en, e1, e2, e3, e4):
    if abs(e1 - e3) > _TOL and abs(e4 - e2) > _TOL: return (e1 - en) * (en - e2) / ((e1 - e3) * (e4 - e2))
    s = _regularize(min(abs(e3 - e1), abs(e4 - e2)))
    num = (e1 - en + s) * (en - e2 + s)
    den = (e1 - e3 + s) * (e4 - e2 + s)
    return float(np.real(num / den))

def _K2(en, e1, e2, e3):
    if abs(e1 - e3) > _TOL and abs(e1 - e2) > _TOL: return (en - e1) / ((e2 - e1) * (e3 - e1))
    s = _regularize(min(abs(e3 - e1), abs(e1 - e2)))
    num = (en - e1 + s)
    den = (e2 - e1 + s) * (e3 - e1 + s)
    return float(np.real(num / den))

def _K1(en, e1, e2):
    if abs(e1 - e2) > _TOL: return (e1 - en) / ((e2 - e1) * (e2 - e1))
    s = _regularize(abs(e1 - e2))
    num = (e1 - en + s)
    den = (e2 - e1 + s) * (e2 - e1 + s)
    return float(np.real(num / den))

def _dos_reorder(en, e):
    # Returns (flag, order, sorted_e)
    order = np.argsort(e)
    se = e[order].copy()

    if (se[0] <= en <= se[3]) and abs(se[3] - se[0]) < _TOL: return 6, order, se
    if se[0] <= en <= se[1]: return 1, order, se
    if se[1] <= en <= se[2]: return 2, order, se
    if se[2] <= en <= se[3]: return 3, order, se
    if en < se[0]: return 4, order, se
    if se[3] < en: return 5, order, se

    return -1, order, se

def _fun_case1(en, e):
    e1, e2, e3, e4 = e
    ci = np.zeros(4, dtype=float)
    ci[0] = _K2(en, e1, e2, e4) * _F(en, e2, e1, e1, e3) \
          + _K2(en, e1, e2, e3) * _F(en, e3, e1, e1, e4) \
          + _K2(en, e1, e3, e4) * _F(en, e4, e1, e1, e2)
    ci[1] = -_K1(en, e1, e2) * _F(en, e1, e1, e3, e4)
    ci[2] = -_K1(en, e1, e3) * _F(en, e1, e1, e2, e4)
    ci[3] = -_K1(en, e1, e4) * _F(en, e1, e1, e2, e3)
    return ci

def _fun_case2(en, e):
    e1, e2, e3, e4 = e
    ci = np.zeros(4, dtype=float)
    ci[0] = 0.5 * (_K1(en, e3, e1) * (
                    _F(en, e3, e2, e2, e4) +
                    _F(en, e4, e1, e2, e4) +
                    _F(en, e3, e1, e2, e4)) +
                   _K1(en, e4, e1) * (
                    _F(en, e4, e1, e2, e3) +
                    _F(en, e4, e2, e2, e3) +
                    _F(en, e3, e1, e2, e3)))
    ci[1] = 0.5 * (_K1(en, e3, e2) * (
                    _F(en, e3, e2, e1, e4) +
                    _F(en, e4, e2, e1, e4) +
                    _F(en, e3, e1, e1, e4)) +
                   _K1(en, e4, e2) * (
                    _F(en, e3, e2, e1, e3) +
                    _F(en, e4, e1, e1, e3) +
                    _F(en, e4, e2, e1, e3)))
    ci[2] = 0.5 * (-_K1(en, e2, e3) * (
                    _F(en, e3, e2, e1, e4) +
                    _F(en, e4, e2, e1, e4) +
                    _F(en, e3, e1, e1, e4)) -
                   _K1(en, e1, e3) * (
                    _F(en, e3, e2, e2, e4) +
                    _F(en, e4, e1, e2, e4) +
                    _F(en, e3, e1, e2, e4)))
    ci[3] = 0.5 * (-_K1(en, e2, e4) * (
                    _F(en, e3, e2, e1, e3) +
                    _F(en, e4, e1, e1, e3) +
                    _F(en, e4, e2, e1, e3)) -
                   _K1(en, e1, e4) * (
                    _F(en, e4, e1, e2, e3) +
                    _F(en, e4, e2, e2, e3) +
                    _F(en, e3, e1, e2, e3)))
    return ci

def _fun_case3(en, e):
    e1, e2, e3, e4 = e
    ci = np.zeros(4, dtype=float)
    ci[0] =  _K1(en, e4, e1) * _F(en, e4, e4, e2, e3)
    ci[1] =  _K1(en, e4, e2) * _F(en, e4, e4, e1, e3)
    ci[2] =  _K1(en, e4, e3) * _F(en, e4, e4, e1, e2)
    ci[3] = -_K2(en, e4, e3, e1) * _F(en, e4, e3, e2, e4) \
            -_K2(en, e4, e2, e3) * _F(en, e4, e2, e1, e4) \
            -_K2(en, e4, e1, e2) * _F(en, e4, e1, e3, e4)
    return ci

def _dos_corner_weights(en, e):
    flag, order, se = _dos_reorder(en, e)
    if   flag == 1: ci = _fun_case1(en, se)
    elif flag == 2: ci = _fun_case2(en, se)
    elif flag == 3: ci = _fun_case3(en, se)
    elif flag in (4, 5): 
        ci = np.zeros(4, dtype=float)
    elif flag == 6:
        ci = np.full(4, 0.25, dtype=float)
    else: raise ValueError("Unexpected flag in tetra reorder")
    return flag, order, ci

def dos_tetra_weights_3d(eigenvalues, energy, k_points):
    """
    Pure-Python version of dos_tetra_weights_3d.
    Inputs:
      - eigenvalues: 1D ndarray, band energies for each k-point (one band)
      - energy: float, evaluation energy
      - k_points: int ndarray with shape (5, ntet); corners are rows 1..4
    Returns:
      - cti: float ndarray (4, ntet), corner weights per tetrahedron
    """
    eigk = np.asarray(eigenvalues, dtype=float)
    itt = np.asarray(k_points, dtype=np.int64)
    if itt.ndim != 2 or itt.shape[0] != 5:
        raise ValueError("k_points must have shape (5, ntet)")
    ntet = itt.shape[1]
    cti = np.zeros((4, ntet), dtype=float)

    for it in range(ntet):
        # rows 1..4 index the four corners
        corners = itt[1:5, it].astype(np.int64)
        e = eigk[corners].astype(float).copy()
        _, order, ci = _dos_corner_weights(energy, e)
        # Map sorted corner weights back to original corner ordering
        # order[j] is original corner index 0..3 for sorted position j
        cti[order, it] = ci
    return cti
