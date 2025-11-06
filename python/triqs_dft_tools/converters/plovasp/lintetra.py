import numpy as np

_TOL = 1e-8

def _regularize_scalar(value):
    """Regularize a scalar to avoid numerical underflow or division by zero."""
    return max(float(value), 1.0e-20) / 100.0

def _select(mask, *arrays):
    """Helper to select masked values for multiple arrays efficiently."""
    idx = np.nonzero(mask)[0]
    return idx, [a[idx] for a in arrays]

# === Auxiliary numerical functions (formerly F, K1, K2) ===

def _stable_fraction_F(e_eval, e1, e2, e3, e4):
    """
    Stable evaluation of the F-term used in tetrahedron DOS integration.

    Parameters
    ----------
    e_eval : float or ndarray
        Evaluation energy.
    e1, e2, e3, e4 : ndarray
        Corner energies of the tetrahedron.

    Returns
    -------
    ndarray
        Evaluated F-term.
    """
    e1, e2, e3, e4, e_eval = map(np.asarray, (e1, e2, e3, e4, e_eval))
    mask = (np.abs(e1 - e3) > _TOL) & (np.abs(e4 - e2) > _TOL)
    safe_val = ((e1 - e_eval) * (e_eval - e2)) / ((e1 - e3) * (e4 - e2))

    s_real = np.maximum(np.minimum(np.abs(e3 - e1), np.abs(e4 - e2)), 1.0e-20) / 100.0
    s = 1j * s_real
    num = (e1 - e_eval + s) * (e_eval - e2 + s)
    den = (e1 - e3 + s) * (e4 - e2 + s)
    fallback = (num / den).real
    return np.where(mask, safe_val, fallback.astype(np.float64))


def _stable_fraction_K1(e_eval, e1, e2):
    """
    Stable evaluation of the first K-term used in tetrahedron DOS integration.

    Parameters
    ----------
    e_eval : float or ndarray
        Evaluation energy.
    e1, e2 : ndarray
        Corner energies of the tetrahedron.

    Returns
    -------
    ndarray
        Evaluated K1-term.
    """
    e1, e2, e_eval = map(np.asarray, (e1, e2, e_eval))
    mask = np.abs(e1 - e2) > _TOL
    safe_val = (e1 - e_eval) / ((e2 - e1) ** 2)

    s_real = np.maximum(np.abs(e1 - e2), 1.0e-20) / 100.0
    s = 1j * s_real
    num = e1 - e_eval + s
    den = (e2 - e1 + s) ** 2
    fallback = (num / den).real
    return np.where(mask, safe_val, fallback.astype(np.float64))


def _stable_fraction_K2(e_eval, e1, e2, e3):
    """
    Stable evaluation of the second K-term used in tetrahedron DOS integration.

    Parameters
    ----------
    e_eval : float or ndarray
        Evaluation energy.
    e1, e2, e3 : ndarray
        Corner energies of the tetrahedron.

    Returns
    -------
    ndarray
        Evaluated K2-term.
    """
    e1, e2, e3, e_eval = map(np.asarray, (e1, e2, e3, e_eval))
    mask = (np.abs(e1 - e3) > _TOL) & (np.abs(e1 - e2) > _TOL)
    safe_val = (e_eval - e1) / ((e2 - e1) * (e3 - e1))

    s_real = np.maximum(np.minimum(np.abs(e3 - e1), np.abs(e1 - e2)), 1.0e-20) / 100.0
    s = 1j * s_real
    num = e_eval - e1 + s
    den = (e2 - e1 + s) * (e3 - e1 + s)
    fallback = (num / den).real
    return np.where(mask, safe_val, fallback.astype(np.float64))


# === Main driver ===

def dos_tetra_weights_3d(eigenvalues, energy, k_points):
    """
    Compute tetrahedron corner weights for 3D DOS integration.

    This version is fully vectorized in NumPy, operating on all tetrahedra
    simultaneously without MPI or Python loops.

    Parameters
    ----------
    eigenvalues : (n_kpoints,) array_like of float
        Energies at k-points.
    energy : float
        Target energy for DOS evaluation.
    k_points : (5, n_tetra) array_like of int
        Tetrahedron connectivity. Only rows [1:5] are used for corner indices.

    Returns
    -------
    corner_weights : (4, n_tetra) ndarray of float
        Corner weights for each tetrahedron at the specified energy.
    """
    eigk = np.asarray(eigenvalues, dtype=float)
    tetra = np.asarray(k_points, dtype=np.int64)
    if tetra.ndim != 2 or tetra.shape[0] != 5:
        raise ValueError("k_points must have shape (5, n_tetra)")

    n_tetra = tetra.shape[1]
    corners = tetra[1:5, :]  # (4, n_tetra)
    corner_energies = eigk[corners]

    # Sort each tetrahedron's corner energies ascending
    order = np.argsort(corner_energies, axis=0)
    sorted_energies = np.take_along_axis(corner_energies, order, axis=0)

    e1, e2, e3, e4 = sorted_energies
    e_eval = float(energy)

    # Determine which energy range each tetrahedron falls into
    flag_uniform = (e1 <= e_eval) & (e_eval <= e4) & (np.abs(e4 - e1) < _TOL)
    flag_case1 = (e1 <= e_eval) & (e_eval <= e2) & (~flag_uniform)
    flag_case2 = (e2 <= e_eval) & (e_eval <= e3)
    flag_case3 = (e3 <= e_eval) & (e_eval <= e4)

    weights_sorted = np.zeros_like(sorted_energies, dtype=float)
    idx = lambda mask: np.nonzero(mask)[0]

    # === Case 6: uniform tetrahedra (degenerate energies)
    if flag_uniform.any():
        weights_sorted[:, idx(flag_uniform)] = 0.25

    # === Case 1
    if flag_case1.any():
        i, (ge1, ge2, ge3, ge4) = _select(flag_case1, e1, e2, e3, e4)
        ee = e_eval

        w0 = (_stable_fraction_K2(ee, ge1, ge2, ge4) * _stable_fraction_F(ee, ge2, ge1, ge1, ge3)
              + _stable_fraction_K2(ee, ge1, ge2, ge3) * _stable_fraction_F(ee, ge3, ge1, ge1, ge4)
              + _stable_fraction_K2(ee, ge1, ge3, ge4) * _stable_fraction_F(ee, ge4, ge1, ge1, ge2))
        w1 = -_stable_fraction_K1(ee, ge1, ge2) * _stable_fraction_F(ee, ge1, ge1, ge3, ge4)
        w2 = -_stable_fraction_K1(ee, ge1, ge3) * _stable_fraction_F(ee, ge1, ge1, ge2, ge4)
        w3 = -_stable_fraction_K1(ee, ge1, ge4) * _stable_fraction_F(ee, ge1, ge1, ge2, ge3)
        weights_sorted[:, i] = np.vstack([w0, w1, w2, w3])

    # === Case 2
    if flag_case2.any():
        i, (ge1, ge2, ge3, ge4) = _select(flag_case2, e1, e2, e3, e4)
        ee = e_eval

        w0 = 0.5 * (
            _stable_fraction_K1(ee, ge3, ge1)
            * (_stable_fraction_F(ee, ge3, ge2, ge2, ge4)
               + _stable_fraction_F(ee, ge4, ge1, ge2, ge4)
               + _stable_fraction_F(ee, ge3, ge1, ge2, ge4))
            + _stable_fraction_K1(ee, ge4, ge1)
            * (_stable_fraction_F(ee, ge4, ge1, ge2, ge3)
               + _stable_fraction_F(ee, ge4, ge2, ge2, ge3)
               + _stable_fraction_F(ee, ge3, ge1, ge2, ge3))
        )

        w1 = 0.5 * (
            _stable_fraction_K1(ee, ge3, ge2)
            * (_stable_fraction_F(ee, ge3, ge2, ge1, ge4)
               + _stable_fraction_F(ee, ge4, ge2, ge1, ge4)
               + _stable_fraction_F(ee, ge3, ge1, ge1, ge4))
            + _stable_fraction_K1(ee, ge4, ge2)
            * (_stable_fraction_F(ee, ge3, ge2, ge1, ge3)
               + _stable_fraction_F(ee, ge4, ge1, ge1, ge3)
               + _stable_fraction_F(ee, ge4, ge2, ge1, ge3))
        )

        w2 = 0.5 * (
            -_stable_fraction_K1(ee, ge2, ge3)
            * (_stable_fraction_F(ee, ge3, ge2, ge1, ge4)
               + _stable_fraction_F(ee, ge4, ge2, ge1, ge4)
               + _stable_fraction_F(ee, ge3, ge1, ge1, ge4))
            - _stable_fraction_K1(ee, ge1, ge3)
            * (_stable_fraction_F(ee, ge3, ge2, ge2, ge4)
               + _stable_fraction_F(ee, ge4, ge1, ge2, ge4)
               + _stable_fraction_F(ee, ge3, ge1, ge2, ge4))
        )

        w3 = 0.5 * (
            -_stable_fraction_K1(ee, ge2, ge4)
            * (_stable_fraction_F(ee, ge3, ge2, ge1, ge3)
               + _stable_fraction_F(ee, ge4, ge1, ge1, ge3)
               + _stable_fraction_F(ee, ge4, ge2, ge1, ge3))
            - _stable_fraction_K1(ee, ge1, ge4)
            * (_stable_fraction_F(ee, ge4, ge1, ge2, ge3)
               + _stable_fraction_F(ee, ge4, ge2, ge2, ge3)
               + _stable_fraction_F(ee, ge3, ge1, ge2, ge3))
        )

        weights_sorted[:, i] = np.vstack([w0, w1, w2, w3])

    # === Case 3
    if flag_case3.any():
        i, (ge1, ge2, ge3, ge4) = _select(flag_case3, e1, e2, e3, e4)
        ee = e_eval

        w0 = _stable_fraction_K1(ee, ge4, ge1) * _stable_fraction_F(ee, ge4, ge4, ge2, ge3)
        w1 = _stable_fraction_K1(ee, ge4, ge2) * _stable_fraction_F(ee, ge4, ge4, ge1, ge3)
        w2 = _stable_fraction_K1(ee, ge4, ge3) * _stable_fraction_F(ee, ge4, ge4, ge1, ge2)
        w3 = (
            -_stable_fraction_K2(ee, ge4, ge3, ge1) * _stable_fraction_F(ee, ge4, ge3, ge2, ge4)
            - _stable_fraction_K2(ee, ge4, ge2, ge3) * _stable_fraction_F(ee, ge4, ge2, ge1, ge4)
            - _stable_fraction_K2(ee, ge4, ge1, ge2) * _stable_fraction_F(ee, ge4, ge1, ge3, ge4)
        )

        weights_sorted[:, i] = np.vstack([w0, w1, w2, w3])

    # === Remap to original corner order
    corner_weights = np.zeros_like(corner_energies, dtype=float)
    corner_weights[order, np.arange(n_tetra)[None, :]] = weights_sorted

    return corner_weights
