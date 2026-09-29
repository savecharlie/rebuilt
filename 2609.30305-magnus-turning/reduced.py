"""The reduced scalar problem, Eq. (8)-(11) of Gaur 2609.30305.

    dQ/du = -lambda e^{2 mu u} / sqrt(1+Q^2),   Q(0) = tan(theta)
    closure   H = int_0^Psi sin gamma du = 0      (returns to launch height)
    chord     = (1/kL) int_0^Psi cos gamma e^{iu} du

Carried as one ODE system so the quadratures share the integrator's error
control with the state they depend on.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq


def _rhs(u, y, lam, mu):
    Q = y[0]
    rt = np.sqrt(1.0 + Q*Q)
    g = Q/rt          # sin gamma
    c = 1.0/rt        # cos gamma
    return [-lam*np.exp(2*mu*u)/rt, g, c*np.cos(u), c*np.sin(u)]


def arc(lam, mu, theta, Psi=np.pi, rtol=1e-13, atol=1e-15):
    s = solve_ivp(_rhs, (0.0, Psi), [np.tan(theta), 0.0, 0.0, 0.0],
                  args=(lam, mu), method="DOP853", rtol=rtol, atol=atol,
                  dense_output=True)
    return s


def closure_lambda(mu, theta, Psi=np.pi, **kw):
    """lambda such that the ball is back at release height after turn Psi."""
    f = lambda lam: arc(lam, mu, theta, Psi, **kw).y[1, -1]
    lo = 1e-9
    hi = max(1.0, 4*np.tan(theta) if theta > 0 else 1.0)
    while f(hi) > 0:
        hi *= 2
        if hi > 1e9:
            raise RuntimeError("no closure bracket")
    return brentq(f, lo, hi, xtol=1e-16, rtol=8.9e-16, maxiter=300)


def chord_bearing_deficit(mu, theta, Psi=np.pi, **kw):
    """pi/2 - arg(chord) in radians, at the closure lambda."""
    lam = closure_lambda(mu, theta, Psi, **kw)
    y = arc(lam, mu, theta, Psi, **kw).y[:, -1]
    return np.pi/2 - np.arctan2(y[3], y[2]), lam
