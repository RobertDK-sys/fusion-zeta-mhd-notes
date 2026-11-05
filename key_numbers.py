import numpy as np, sys
sys.path.insert(0,'.')
from icrf_core import *
from loading import rc

CASES = {
 'SPARC': dict(B=12.2, R0=1.85, a=0.57, ne=4.0e20, f=120e6, kpar=15.0,
               fr3=0.04, nsep=5e19, w=0.05, Lp=0.25, Z0=20.0),
 'ITER':  dict(B=5.3,  R0=6.2,  a=2.0,  ne=1.0e20, f=53e6,  kpar=8.0,
               fr3=0.04, nsep=4e19, w=0.15, Lp=0.30, Z0=20.0),
}
comp = lambda ne, fr3: dict(D=(1-fr3)/2*ne, T=(1-fr3)/2*ne, He3=fr3*ne, e=ne*(1+fr3))

print("=== (1) Onde rapide : vitesse d'Alfven et nombre d'onde perpendiculaire ===")
for tag, c in CASES.items():
    fr = {'D':(1-c['fr3'])/2, 'T':(1-c['fr3'])/2, 'He3':c['fr3']}
    mbar = sum(fr[k]*mass(k) for k in fr)
    va = v_alfven(c['B'], c['ne'], mbar); w = 2*np.pi*c['f']
    n2 = nperp2(w, c['B'], comp(c['ne'], c['fr3']), 0.0).real
    kperp = np.sqrt(n2)*w/299792458.0
    rho = np.sqrt(2*20000*1.602176634e-19/mass('D'))/(2*np.pi*fcyc('D', c['B']))
    print(f"  {tag}: v_A={va:.3e} m/s  c/v_A={299792458.0/va:.1f}  w/v_A={w/va:.1f} /m"
          f"  k_perp(cold,exact)={kperp:.1f} /m  k_perp*a={kperp*c['a']:.1f}"
          f"  ratio={kperp/(w/va):.2f}  k_perp*rho_D(T=20keV)={kperp*rho:.2f}")

print("\n=== (2) Coupure de densite de l'onde rapide (S=n_par^2) ===")
for tag, c in CASES.items():
    w = 2*np.pi*c['f']; npar = c['kpar']*299792458.0/w
    ne = np.logspace(15, np.log10(c['ne']*2), 4000)
    n2 = nperp2(w, c['B'], comp(ne, c['fr3']), npar).real
    idx = np.where(np.sign(n2[1:]) != np.sign(n2[:-1]))[0]
    ncut = ne[idx[0]+1] if len(idx) else np.nan
    print(f"  {tag}: n_par={npar:.2f} (k_par={c['kpar']:.0f}/m)  n_cut={ncut:.2e} = {ncut/c['ne']:.1%} de n_e0")
    for kk in [6.0, 10.0, 12.0]:
        npk = kk*299792458.0/w
        n2k = nperp2(w, c['B'], comp(ne, c['fr3']), npk).real
        ik = np.where(np.sign(n2k[1:]) != np.sign(n2k[:-1]))[0]
        print(f"      k_par={kk:.0f}/m -> n_cut={ne[ik[0]+1]:.2e} ({ne[ik[0]+1]/c['ne']:.1%})" if len(ik) else "")

print("\n=== (3) Couche de resonance et modulation en frequence ===")
for tag, c in CASES.items():
    fc = fcyc('He3', c['B']); Rres = c['R0']*fc/c['f']
    print(f"  {tag}: f_c(3He, B0)={fc/1e6:.2f} MHz ; a {c['f']/1e6:.0f} MHz: R_res={Rres:.3f} m"
          f" (axe {c['R0']} m, soit {Rres-c['R0']:+.3f} m)")
    for df in [-1e6, 1e6]:
        print(f"      f{df/1e6:+.0f} MHz -> R_res={c['R0']*fc/(c['f']+df):.3f} m (decalage {(c['R0']*fc/(c['f']+df)-Rres)*100:+.2f} cm)")
    vt = np.sqrt(2*20000*1.602176634e-19/mass('He3')); w = 2*np.pi*c['f']
    kpar = c['kpar']; dxdop = kpar*vt*Rres/w
    print(f"      largeur Doppler k_par*v_t/|dw/dx| = {dxdop*100:.1f} cm   (v_t(3He,20keV)={vt:.2e} m/s)")

print("\n=== (4) R_c(d) et loi exponentielle (modele 1-D, geometrie declaree) ===")
for tag, c in CASES.items():
    fr = {'D':(1-c['fr3'])/2,'T':(1-c['fr3'])/2,'He3':c['fr3']}
    ds = np.arange(0.02, 0.31, 0.01); out=[]
    for d in ds:
        out.append(rc(d, fr, c['B'], c['f'], c['kpar'], c['w'], c['Lp'], c['nsep'], 0.02, c['ne'])[0])
    out = np.array(out); ok = out>0
    p = np.polyfit(ds[ok], np.log(out[ok]), 1)
    print(f"  {tag}: R_c(d=14cm)={np.exp(p[1]+p[0]*0.14):.2f} ohm, R_c(d=6cm)={np.exp(p[1]+p[0]*0.06):.2f}, "
          f"pente={-p[0]:.1f} /m (2k_par={2*c['kpar']:.0f}) ; R_c(30cm)/R_c(6cm)={np.exp(p[0]*0.24):.2e}")

print("\n=== (5) Limite de tension : P_max = V^2 R_c /(2 Z0^2) par ligne ===")
for V in [45e3, 60e3]:
    row = [f"Z0={z:.0f}: P_max={V**2*1.0/(2*z*z)/1e6:.2f} MW a R_c=1 ohm" for z in [15.,20.,30.]]
    print(f"  V_max={V/1e3:.0f} kV -> " + " ; ".join(row))

print("\n=== (6) Coupures ECRH (O-mode) ===")
for f in [140e9, 170e9, 240e9]:
    print(f"  f={f/1e9:.0f} GHz : n_cut = 1.24e-2 f^2 = {1.24e-2*f**2:.2e} m^-3")
print(f"  pour n_e=4e20 m^-3 il faut f >= {np.sqrt(4e20/1.24e-2)/1e9:.0f} GHz (SPARC) ; "
      f"n_e=1e20 -> {np.sqrt(1e20/1.24e-2)/1e9:.0f} GHz (ITER)")
