"""Fast-wave (ICRF) core quantities: dispersion, cutoffs, antenna loading.

Conventions: time dependence exp(-i w t); B0 along z (toroidal);
propagation along x (radial); k_parallel = n_par * w / c along z.
Cold plasma Stix parameters R, L, S, D, P; fast-wave local dispersion
    n_perp^2 = (R - n_par^2)(L - n_par^2) / (S - n_par^2).
"""
import numpy as np

e, me, mp, eps0, mu0, c = 1.602176634e-19, 9.1093837015e-31, 1.67262192369e-27, 8.8541878128e-12, 1.25663706212e-6, 299792458.0

IONS = {  # name: (charge state, mass in u)
    'e': (-1, 5.48579909e-4),
    'H': (1, 1.007276), 'D': (1, 2.013553), 'T': (1, 3.015500),
    'He3': (2, 3.014932), 'He4': (2, 4.001506),
}

def mass(name):
    return IONS[name][1] * 1.66053906660e-27

def fcyc(name, B):
    """cyclotron frequency in Hz"""
    Z, A = IONS[name]
    return Z * e * B / (2 * np.pi * mass(name))

def wpl2(name, n):
    """squared plasma frequency (rad/s)^2 for density n (m^-3)"""
    Z, A = IONS[name]
    return n * Z * Z * e * e / (eps0 * mass(name))

def stix(w, B, comp):
    """comp: dict species -> density (m^-3). Returns R, L, S, D, P."""
    w, n0 = np.broadcast_arrays(np.asarray(w, float), np.asarray(list(comp.values())[0], float))
    shape = w.shape
    R = np.ones(shape, dtype=complex)
    L = np.ones(shape, dtype=complex)
    P = np.ones(shape, dtype=complex)
    for sp, n in comp.items():
        Z, A = IONS[sp]
        O = 2 * np.pi * fcyc(sp, B)
        wp2 = np.broadcast_to(wpl2(sp, np.asarray(n, float)), shape)
        R -= wp2 / (w * (w + O))
        L -= wp2 / (w * (w - O))
        P -= wp2 / (w * w)
    S = 0.5 * (R + L); D = 0.5 * (R - L)
    return R, L, S, D, P

def nperp2(w, B, comp, npar):
    R, L, S, D, P = stix(w, B, comp)
    return (R - npar**2) * (L - npar**2) / (S - npar**2)

def v_alfven(B, n, mi):
    return B / np.sqrt(mu0 * n * mi)

def ncut_fast(w, B, comp_template, npar, lo=1e13, hi=1e22):
    """density scale s such that S = npar^2 (fast-wave cutoff), by bisection on scale factor."""
    f = lambda s: stix(w, B, {k: v * s for k, v in comp_template.items()})[2].real - npar**2
    a, b = lo, hi
    if f(a) > 0:
        return None            # already propagating at 1e13: cutoff below this
    if f(b) < 0:
        return None            # never reaches cutoff below 1e22
    for _ in range(200):
        m = np.sqrt(a * b)
        if f(m) > 0: b = m
        else: a = m
    return np.sqrt(a * b)
