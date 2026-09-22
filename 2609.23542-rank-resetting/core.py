"""
Independent rebuild of arXiv:2609.23542 (Vatash, Goldfarb, Rudyak, Roichman)
"Rank-dependent optimal resetting in multiparticle search"

Two instruments, deliberately sharing nothing:

  A. DETERMINISTIC.  Numerical Laplace inversion of the renewal relation
     (their Eq. 1-2) -> Q_r(t), then quadrature for <T_(k)> via Eq. 3-4.

  B. EXACT MONTE CARLO.  No time discretisation at all.  A single searcher's
     free first-passage time is Levy(c) with c = d^2/(2D), sampled exactly as
     c/Z^2, Z~N(0,1); resetting is an exact renewal loop over Exp(r) clocks.

Units throughout: D = 1.  Distances in units of d, rates in units of D/d^2.

Iris, Sep 22 2026.
"""
import numpy as np
from math import erf, exp, sqrt, pi, log
from scipy.special import erf as serf
from scipy.optimize import brentq, minimize_scalar
from scipy.integrate import quad
from scipy.special import comb

D = 1.0

# ---------------------------------------------------------------- free search
def Q0(t, d=1.0):
    """Survival of a free 1D Brownian searcher started distance d from the
    absorbing target.  Q0 = erf(d / sqrt(4 D t)).  (Laplace: (1-e^{-d sqrt(s/D)})/s)"""
    t = np.asarray(t, dtype=float)
    out = np.ones_like(t)
    m = t > 0
    out[m] = serf(d / np.sqrt(4.0 * D * t[m]))
    return out if out.ndim else float(out)

def Q0t_L(s, d=1.0):
    """Laplace transform of Q0, their Eq. 2."""
    return (1.0 - np.exp(-d * np.sqrt(s / D))) / s

# ------------------------------------------------- A: deterministic instrument
def Qr_L(s, r, d=1.0):
    """Laplace transform of the survival WITH resetting, their Eq. 1."""
    q = Q0t_L(s + r, d)
    return q / (1.0 - r * q)

def mean_fpt_1p(r, d=1.0):
    """Single-particle mean FPT under resetting = Qr_L(0).  Closed form:
    (1/r)(exp(d sqrt(r/D)) - 1).  Evans & Majumdar 2011."""
    if r <= 0:
        return np.inf
    return (np.exp(d * np.sqrt(r / D)) - 1.0) / r

def Qr_time(ts, r, d=1.0, degree=None):
    """Q_r(t) by numerical Laplace inversion (Talbot) of Eq. 1.  mpmath."""
    import mpmath as mp
    old = mp.mp.dps
    mp.mp.dps = 30
    dd = mp.mpf(d)
    def F(s):
        # Eq.1-2 combined and simplified.  NOTE: Q0 is evaluated at s+r, not s.
        # Writing sqrt(s/D) here instead of sqrt((s+r)/D) is a typo that still
        # returns finite plausible numbers -- it cost an hour on Sep 22 2026.
        w = mp.sqrt((s + r) / D)
        return (1 - mp.e**(-dd * w)) / (s + r * mp.e**(-dd * w))
    out = []
    for t in np.atleast_1d(ts):
        out.append(float(mp.invertlaplace(F, mp.mpf(float(t)), method='talbot',
                                          degree=degree or 24)))
    mp.mp.dps = old
    return np.array(out)

# ------------------------------------------------------- B: exact MC instrument
def sample_fpt(r, d=1.0, size=1, rng=None):
    """EXACT samples of the first-passage time of one resetting searcher.
    No dt.  Free FPT ~ Levy(c), c = d^2/(2D), sampled as c/Z^2.
    Renewal loop: draw a reset clock R~Exp(r) and a free attempt tau;
    if tau < R the search ends, else the elapsed time is R and we restart."""
    rng = rng or np.random.default_rng()
    c = d * d / (2.0 * D)
    n = int(size)
    elapsed = np.zeros(n)
    alive = np.arange(n)
    while alive.size:
        m = alive.size
        z = rng.standard_normal(m)
        tau = c / (z * z)
        if r > 0:
            R = rng.exponential(1.0 / r, m)
        else:
            R = np.full(m, np.inf)
        done = tau < R
        elapsed[alive[done]] += tau[done]
        elapsed[alive[~done]] += R[~done]
        alive = alive[~done]
    return elapsed

def mc_ordered(r, ds, ntrial, rng=None):
    """MC estimate of <T_(k)> for k=1..N.  ds = list of initial distances."""
    rng = rng or np.random.default_rng()
    N = len(ds)
    T = np.empty((ntrial, N))
    for i, d in enumerate(ds):
        T[:, i] = sample_fpt(r, d=d, size=ntrial, rng=rng)
    T.sort(axis=1)
    return T.mean(axis=0), T.std(axis=0, ddof=1) / sqrt(ntrial)
