"""
How BIG does the imprint of fitness get?

The theorem says the imprint vanishes and decays at the spectral gap.  It does
not say how large it gets first.  Claim (mine):

In the no-mixing limit the composition p = u/U obeys dp/ds = (r - rbar)p with
ds = g(U)dt, so p(s) = psi e^{rs}/M(s) with M(s) = \int psi e^{rs} dx, and
d ln U/ds = rbar = K'(s), K = ln M.  Integrating, ln(U/U0) = K(s).  Growth stops
at U*, so the TOTAL selection time is

        s*  solves  K(s*) = ln(U*/U0),      K = cumulant generating function
                                                 of fitness r under psi

and the largest imprint any switching dynamics can fail to prevent is

        p_ceiling(x) = psi(x) e^{r(x)s*} / M(s*)        (exponential tilt of psi)

Note s* has the SIGN of ln(U*/U0): a population that shrinks to carrying
capacity tilts towards LOW fitness.
"""
import numpy as np
from scipy.optimize import brentq


def selection_time(psi, r, h, growth_log):
    """s* solving K(s*) = growth_log, K(s)=ln \int psi e^{rs}."""
    def K(s):
        m = (r*s).max()
        return m + np.log(h*np.sum(psi*np.exp(r*s - m)))
    lo, hi = -1e3, 1e3
    # bracket
    a, b = -1.0, 1.0
    while K(b) < growth_log and b < 1e4: b *= 2
    while K(a) > growth_log and a > -1e4: a *= 2
    return brentq(lambda s: K(s) - growth_log, a, b, xtol=1e-14, rtol=1e-15)


def ceiling_composition(psi, r, h, growth_log):
    s = selection_time(psi, r, h, growth_log)
    m = (r*s).max()
    w = psi*np.exp(r*s - m)
    return w/(h*np.sum(w)), s
