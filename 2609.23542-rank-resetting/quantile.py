"""
Third, independent road to phi_c: no perturbation theory at all.
Compute the phi-quantile t_phi(r) of the restarted first-passage distribution
directly from the Laplace inversion, and look for an interior minimum.
Prediction from the closed form: interior minimum iff phi > 0.4123101755.
"""
import numpy as np, core
from scipy.optimize import brentq

def t_quantile(phi, r, d=1.0):
    """t such that F_r(t) = phi, i.e. Q_r(t) = 1-phi."""
    f = lambda lt: core.Qr_time([np.exp(lt)], r, d=d, degree=30)[0] - (1-phi)
    lo, hi = -8.0, 8.0
    while f(hi) > 0: hi += 3
    return np.exp(brentq(f, lo, hi, xtol=1e-10))

if __name__ == "__main__":
    rs = np.array([1e-4, 0.02, 0.05, 0.1, 0.2, 0.4, 0.8, 1.6])
    print("t_phi(r) / t_phi(0) -- interior minimum means resetting helps that quantile")
    print(f"{'phi':>7} " + " ".join(f"{r:>9.4g}" for r in rs) + "   min at")
    for phi in [0.30, 0.38, 0.405, 0.4123101755, 0.42, 0.50, 0.70]:
        vals = np.array([t_quantile(phi, r) for r in rs])
        rel = vals/vals[0]
        i = int(np.argmin(rel))
        tag = "r=0 (none)" if i == 0 else f"r~{rs[i]:.3g}  INTERIOR"
        print(f"{phi:7.4f} " + " ".join(f"{v:9.6f}" for v in rel) + f"   {tag}")
    print()
    print("predicted threshold phi_c = 0.4123101755")
