"""Rebuild PAMGuard's 3D hyperbolic solve (HyperbolicLocaliser.processTOADs3D,
Gillette & Silverman 2008 form) in the geometry of arXiv:2609.20350 and ask:
does sub-sample TDOA noise alone produce km-scale errors?  Iris, fire 282."""
import numpy as np
from scipy.optimize import least_squares

G = np.array([[0,0,-1000],[500,500,-900],[-500,500,-850],[-500,-500,-800],[500,-500,-970]], float)
C = 1450.0

def pamguard_gs(g, t, c=C, weight=True):
    """Line-for-line port of processTOADs3D. t = arrival times; delays[i][j] = t_j - t_i."""
    n = len(g); cen = g.mean(0)
    L, R = [], []
    for i in range(n):
        for j in range(i+1, n):
            delay = -(t[j]-t[i]) * c
            row = np.zeros(3+n-1); rhs = delay**2
            row[:3] = g[i]-g[j]
            rhs += np.sum((g[i]-cen)**2) - np.sum((g[j]-cen)**2)
            row[3+i] = delay
            rhs *= 0.5
            w = np.linalg.norm(g[i]-g[j]) if weight else 1.0   # w = 1/(1/sqrt(scale)) as written
            L.append(row*w); R.append(rhs*w)
    sol, *_ = np.linalg.lstsq(np.array(L), np.array(R), rcond=None)
    return cen + sol[:3], np.array(L)

def nls(g, t, c=C):
    """Plain TDOA least squares (what a simplex/ML localiser converges to)."""
    d = (t[1:]-t[0])*c
    f = lambda s: (np.linalg.norm(g[1:]-s,axis=1)-np.linalg.norm(g[0]-s)) - d
    best = None
    for s0 in [g.mean(0), g.mean(0)+[0,0,500], g.mean(0)+[300,-300,300]]:
        r = least_squares(f, s0)
        if best is None or r.cost < best.cost: best = r
    return best.x

rng = np.random.default_rng(282)
W = np.column_stack([rng.uniform(-500,500,100), rng.uniform(-500,500,100), rng.uniform(-1000,0,100)])

def run(sigma, weight=True):
    eg, en = [], []
    for s in W:
        t = np.linalg.norm(G-s,axis=1)/C + rng.normal(0, sigma, 5)
        eg.append(np.linalg.norm(pamguard_gs(G,t,weight=weight)[0]-s))
        en.append(np.linalg.norm(nls(G,t)-s))
    return np.array(eg), np.array(en)

if __name__ == "__main__":
    e,_ = run(0.0)
    print(f"INSTRUMENT CHECK sigma=0: max GS error {e.max():.2e} m (must be ~1e-9)")
    _, L = pamguard_gs(G, np.linalg.norm(G-W[0],axis=1)/C)
    print(f"cond(L) for whale 0: {np.linalg.cond(L):.2e}")
    for sig in [1e-7, 1e-6, 5e-6, 2e-5]:
        eg, en = run(sig)
        print(f"sigma={sig:.0e}s ({sig*C*100:.2f} cm): GS median {np.median(eg):9.2f}  max {eg.max():10.1f} m | NLS median {np.median(en):7.3f} max {en.max():8.2f} m")
