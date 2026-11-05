"""ICRF antenna loading resistance from a 1-D cold-plasma slab.

B0 || z, strap current sheet K_y on the wall at x=0, k_parallel = n_par k0 along z.
Local dispersion kept exactly:  K(x) = k0^2 (R-n^2)(L-n^2)/(S-n^2).

Wave impedance Z = E_y / (dE_y/dx) obeys the Riccati equation
    dZ/dx = 1 + K(x) Z^2 ,
integrated inward from the plasma core where the wave is outgoing (absorbed),
so no artificial absorbing layer or reflection is introduced.
Antenna boundary condition (B_z jump across the current sheet = mu0 K):
    dE/dx(0+) = i w mu0 K   =>   R_c = -w mu0 Im Z(0) / L_pol
and P = 1/2 R_c I^2 with I = K L_pol.
"""
import numpy as np
from icrf_core import stix, mass, mu0, c, IONS

def profile(x, d, n_sep, lam_sol, n_core, l_core):
    x = np.asarray(x, float); out = np.empty_like(x)
    sol = x < d
    out[sol] = n_sep*np.exp((x[sol]-d)/lam_sol)
    t = np.clip((x[~sol]-d)/l_core, 0, None)
    out[~sol] = n_sep + (n_core-n_sep)*(1-np.exp(-t))
    return out

def K_of(comp_frac, B, w, npar, ne):
    comp = {k: v*ne for k, v in comp_frac.items()}
    comp['e'] = ne*sum(IONS[k][0]*comp_frac[k] for k in comp_frac)
    R, L, S, D, P = stix(w, B, comp)
    den = S - npar**2
    den = np.where(np.abs(den) < 1e-2, den + 1e-2j, den)
    return (w/c)**2*(R-npar**2)*(L-npar**2)/den

def rc(d, comp_frac, B, f, k_par, w_strap, L_pol, n_sep, lam_sol, n_core,
       l_core=0.12, Lx=None, npts=3000):
    w = 2*np.pi*f; k0 = w/c; npar = k_par*c/w
    if Lx is None:
        Lx = d + l_core*(1 + 2.0)      # end well inside the core
    x = np.linspace(0, Lx, npts)
    ne = profile(x, d, n_sep, lam_sol, n_core, l_core)
    K = K_of(comp_frac, B, w, npar, ne)
    Kl = K[-1]
    Z = -1j/np.sqrt(Kl) if Kl.real > 0 else -1.0/np.sqrt(abs(Kl.real))
    dx = x[1]-x[0]
    def f_(Z, xa, Ka):
        return 1.0 + Ka*Z*Z
    for i in range(len(x)-2, -1, -1):
        xa, Kb = x[i], K[i]
        k1 = f_(Z, xa, Kb); k2 = f_(Z+0.5*dx*k1, xa, Kb)
        Zn = Z + dx*k2
        Z = Zn if np.isfinite(Zn) and abs(Zn) < 1e12 else (Z + dx*f_(Z, xa, Kb))
    Rc = -w*mu0*Z.imag*(w_strap/L_pol)   # impedance x strap aspect ratio
    # first fast-wave cut-off from the wall (numerator of the dispersion)
    comp = {k: v*ne for k, v in comp_frac.items()}; comp['e'] = ne
    R, L, S, D, P = stix(w, B, comp)
    num = ((R-npar**2)*(L-npar**2)).real
    idx = np.where(np.sign(num[1:]) != np.sign(num[:-1]))[0]
    xcut = x[idx[0]+1] if len(idx) else np.nan
    return Rc, xcut
