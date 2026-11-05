import numpy as np, sys
sys.path.insert(0,'.')
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from icrf_core import *
from loading import rc

C = {'SPARC': dict(B=12.2, R0=1.85, a=0.57, ne=4.0e20, f=120e6, kpar=15.0,
                   fr3=0.04, nsep=5e19, w=0.05, Lp=0.25),
     'ITER':  dict(B=5.3,  R0=6.2,  a=2.0,  ne=1.0e20, f=53e6,  kpar=8.0,
                   fr3=0.04, nsep=4e19, w=0.15, Lp=0.30)}
comp = lambda ne, fr3: dict(D=(1-fr3)/2*ne, T=(1-fr3)/2*ne, He3=fr3*ne, e=ne*(1+fr3))
fr = lambda fr3: {'D':(1-fr3)/2,'T':(1-fr3)/2,'He3':fr3}

fig = plt.figure(figsize=(13.4, 3.8))
ax0, ax1, ax2 = fig.add_subplot(131), fig.add_subplot(132), fig.add_subplot(133)

# (a) Alfven limit
B, n = 5.3, 1e20; va = v_alfven(B, n, mass('D'))
wr = np.logspace(-2.5, 0.35, 300)*2*np.pi*fcyc('D', B)
ratio = [np.sqrt(max(nperp2(w, B, {'D':n,'e':n}, 0.0).real,1e-30))*va/299792458.0 for w in wr]
ax0.loglog(wr/(2*np.pi*fcyc('D',B)), ratio, lw=1.8)
ax0.axhline(1, color='k', ls='--', lw=0.8)
ax0.axvline(1, color='r', ls=':', lw=0.9)
ax0.annotate('ω = Ω$_i$', (1.2, 1e-2), fontsize=7, color='r')
ax0.set_xlabel(r'$\omega/\Omega_D$'); ax0.set_ylabel(r'$k_\perp/(\omega/v_A)$')
ax0.set_title('(a) $k_\\perp=\\omega/v_A$ exact si $\\omega\\ll\\Omega_i$', fontsize=9)
ax0.grid(alpha=.3, which='both')

# (b) fast-wave index
for tag, c in C.items():
    w = 2*np.pi*c['f']; npar = c['kpar']*299792458.0/w
    ne = np.logspace(15, np.log10(c['ne']*3), 1500)
    cc = comp(ne, c['fr3']); R_, L_, S_, D_, P_ = stix(w, c['B'], cc)
    n2 = nperp2(w, c['B'], cc, npar).real
    y = np.clip(np.sign(n2)*np.sqrt(np.abs(n2)), -25, 60)
    res = np.abs(S_.real - npar**2) < 0.02*npar**2
    y = np.where(res, np.nan, y)
    ax1.semilogx(ne, y, lw=1.6, label=f"{tag} (n$_\\parallel$={npar:.1f})")
    ax1.axvline(c['ne'], color='gray', ls=':', lw=0.7)
ax1.axhline(0, color='k', lw=0.7)
ax1.set_ylim(-25, 60)
ax1.set_xlabel('n$_e$ [m$^{-3}$]'); ax1.set_ylabel(r'sign$(n_\perp^2)\sqrt{|n_\perp^2|}$')
ax1.set_title('(b) évanescente à basse densité,\npropagative en cœur (pointillés = n$_{e0}$)', fontsize=9)
ax1.legend(fontsize=7, loc='upper left'); ax1.grid(alpha=.3)
ax1.annotate('résonance\nhybride ion-ion', (6e18, -22), fontsize=7, color='C1')

# (c) R_c(d)
for tag, c in C.items():
    ds = np.arange(0.02, 0.31, 0.01)
    R = np.array([rc(d, fr(c['fr3']), c['B'], c['f'], c['kpar'], c['w'], c['Lp'],
                      c['nsep'], 0.02, c['ne'])[0] for d in ds])
    ok = R > 0
    p = np.polyfit(ds[ok], np.log(R[ok]), 1)
    ax2.semilogy(ds*100, R, 'o', ms=3, color='C0' if tag=='SPARC' else 'C1')
    ax2.semilogy(ds*100, np.exp(np.polyval(p, ds)), '-',
                 color='C0' if tag=='SPARC' else 'C1',
                 label=f"{tag}: pente {-p[0]:.0f} m$^{{-1}}$ (2k$_\\parallel$={2*c['kpar']:.0f})")
ax2.set_xlabel('distance antenne → couche évanescente  d [cm]')
ax2.set_ylabel('R$_c$ [Ω] (par strap)')
ax2.set_title("(c) R$_c(d)$ : $e^{-2k_\\parallel d/2}$", fontsize=9)
ax2.legend(fontsize=7); ax2.grid(alpha=.3, which='both')
plt.tight_layout(); plt.savefig('results/fig1_fastwave_loading.png', dpi=115)
print('fig1 ok')

fig, ax = plt.subplots(1, 2, figsize=(9.4, 3.7))
Rcv = np.logspace(-2, 1, 300)
for z in [15., 20., 30.]:
    ax[0].loglog(Rcv, (45e3)**2*Rcv/(2*z*z)/1e6, label=f'Z$_0$={z:.0f} Ω')
ax[0].axhline(2.5, color='k', ls='--', lw=0.9)
ax[0].annotate('2,5 MW par ligne\n(20 MW / 8 lignes)', (1.15e-2, 2.7), fontsize=7)
ax[0].set_xlabel('R$_c$ [Ω]'); ax[0].set_ylabel('P$_{max}$ à 45 kV [MW]')
ax[0].set_title(r'$P_{max}=V_{max}^2R_c/2Z_0^2$ : ITER limité par la tension', fontsize=9)
ax[0].legend(fontsize=7); ax[0].grid(alpha=.3, which='both')

f = np.logspace(8, 12, 400)
ax[1].loglog(f/1e9, 1.24e-2*f**2, lw=1.8, label='coupure O-mode $n_c=1.24\\times10^{-2}f^2$')
for lab, ff, nn, mk in [('SPARC 4·10$^{20}$',120,4e20,'s'), ('ITER 10$^{20}$',53,1e20,'o')]:
    ax[1].plot(ff/1e3, nn, mk, color='k'); ax[1].annotate(lab, (ff/1e3, nn),
              textcoords='offset points', xytext=(6,-3), fontsize=7)
for lab, ff in [('140 GHz',140), ('170 GHz',170)]:
    ax[1].axvline(ff/1e3, color='r', ls='--', lw=0.8)
    ax[1].annotate(lab, (ff/1e3, 2e15), fontsize=7, rotation=90, color='r')
ax[1].set_xlabel('f [GHz]'); ax[1].set_ylabel('densité de coupure [m$^{-3}$]')
ax[1].set_title("ECRH : coupure en densité ; l'ICRH n'en a pas", fontsize=9)
ax[1].set_ylim(1e15, 1e23); ax[1].grid(alpha=.3, which='both')
plt.tight_layout(); plt.savefig('results/fig2_voltage_ecrh.png', dpi=115)
print('fig2 ok')
