import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from scipy.special import erf
from closedform import g_closed
from bridge import J

PHI = 0.412310175459; TC = 0.743904342807; NC = 7.32647733; CC = 2.07637

def kc_of_N(N):
    # J rises with k (later ranks benefit more).  A boundary exists only when
    # the FASTEST arrival is already past it, J(N,1) < 0 -- i.e. N > N_c.
    # (First version had this test inverted and drew an empty curve.)
    if J(N, 1.0+1e-9) > 0: return None
    return brentq(lambda k: J(N, k), 1.0+1e-9, N-3.5, xtol=1e-9)

fig, ax = plt.subplots(1, 2, figsize=(11.6, 4.5))

# ---- (a) g(t): the one function both thresholds come from
t = np.geomspace(0.02, 60, 900); g = np.array([g_closed(x) for x in t])
ax[0].axhline(0, color='0.75', lw=.8)
ax[0].plot(t, g, color='#1b3b6f', lw=2)
ax[0].plot([TC], [0], 'o', ms=8, mfc='#d94801', mec='k', zorder=5)
ax[0].set_xscale('log'); ax[0].set_xlabel(r'$t\ \ [d^2/D]$')
ax[0].set_ylabel(r'$g(t)=t\,Q_0(t)-(Q_0*Q_0)(t)$')
ax[0].set_title(r'(a)  $\partial F_r/\partial r|_{r=0}$ — resetting helps where $g>0$', fontsize=10)
ax[0].set_ylim(-0.06, 0.32); ax[0].set_xlim(0.02, 60)
ax[0].annotate(rf'$\tau_c={TC:.6f}$' '\n' rf'$\phi_c=1-Q_0(\tau_c)={PHI:.6f}$',
               xy=(TC, 0), xytext=(2.2, -0.035), fontsize=9,
               arrowprops=dict(arrowstyle='->', color='0.35'))
ax[0].text(0.035, -0.045, 'resetting delays\nthis quantile', fontsize=8, color='0.35')

# ---- (b) the boundary: both published results, one curve
Ns = np.unique(np.round(np.geomspace(8, 4000, 44)).astype(int))
pts = [(N, kc_of_N(N)) for N in Ns]
pts = [(N, k) for N, k in pts if k is not None]
NN = np.array([p[0] for p in pts], float); PH = np.array([p[1] for p in pts])/NN
ax[1].axhline(PHI, color='#d94801', ls='--', lw=1.3)
Nf = np.geomspace(8, 4000, 300)
ax[1].plot(Nf, PHI - CC/Nf, color='0.6', lw=1.2, ls=':')
ax[1].plot(NN, PH, 'o-', color='#1b3b6f', lw=1.6, ms=4.2)
ax[1].plot([NC], [1.0/NC], '*', ms=17, mfc='#2c7b2c', mec='k', zorder=6)
ax[1].set_xscale('log'); ax[1].set_xlabel(r'$N$'); ax[1].set_ylabel(r'$\phi=k/N$')
ax[1].set_title(r'(b)  $J(N,k)=0$ — where a finite $r^*$ appears', fontsize=10)
ax[1].set_xlim(6, 5000); ax[1].set_ylim(0, 0.48)
ax[1].text(200, 0.425, rf'Belan 2020:  $\phi_c={PHI:.6f}$', color='#d94801', fontsize=9)
ax[1].text(150, 0.30, rf'$\phi_c-{CC:.4f}/N$', color='0.45', fontsize=9)
ax[1].text(8.6, 0.055, 'Biroli, Majumdar & Schehr 2023\n' rf'$N_c={NC:.6f}$  ($k=1$)',
           color='#2c7b2c', fontsize=9)
ax[1].text(40, 0.17, 'finite $r^*$ exists\n(above the curve)', fontsize=9, color='#1b3b6f')
ax[1].text(600, 0.06, 'no finite $r^*$', fontsize=9, color='0.4')
plt.tight_layout(); plt.savefig('rank_resetting_rebuild.png', dpi=155)
print('wrote rank_resetting_rebuild.png')
print('boundary sample:', [(int(n), round(p,5)) for n,p in list(zip(NN,PH))[:4]],
      '...', [(int(n), round(p,5)) for n,p in list(zip(NN,PH))[-2:]])
