"""<T_(k)> for N independent resetting searchers, from their Eq.3-4."""
import numpy as np, core
from scipy.special import comb
from scipy.integrate import quad
from scipy.optimize import minimize_scalar

def surv_ordered_identical(Q, N, k):
    """Pr[T_(k)>t] for i.i.d. searchers: fewer than k have arrived."""
    F = 1.0 - Q
    return sum(comb(N, m, exact=True) * F**m * Q**(N-m) for m in range(k))

def mean_Tk(r, N, k, d=1.0, degree=30):
    """<T_(k)> by quadrature over log t, Q_r from Laplace inversion."""
    def f(u):
        t = np.exp(u)
        Q = core.Qr_time([t], r, d=d, degree=degree)[0]
        return surv_ordered_identical(Q, N, k) * t
    hi = np.log(400.0 * core.mean_fpt_1p(r) + 50.0)
    v, _ = quad(f, -18, hi, limit=300)
    return v

def opt_rate(N, k, d=1.0, lo=-6.0, hi=3.5):
    res = minimize_scalar(lambda lr: mean_Tk(np.exp(lr), N, k, d),
                          bounds=(lo, hi), method='bounded',
                          options={'xatol': 1e-5})
    return np.exp(res.x), res.fun
