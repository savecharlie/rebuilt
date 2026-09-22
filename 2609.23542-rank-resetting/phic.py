"""
phi_c from scratch -- no fitting, no simulation, no Laplace inversion.

Large-N quantile limit: <T_(k)> -> t_phi, the phi-quantile of the single-searcher
first-passage distribution, F_r(t_phi) = phi.  A finite optimal reset rate exists
iff t_phi initially DECREASES with r, i.e. d t_phi/dr < 0 at r=0.

Implicit differentiation of F_r(t_phi(r)) = phi:
    d t_phi/dr = -(dF_r/dr) / f(t_phi)          [f > 0]
so the condition is  dF_r/dr|_{r=0} > 0  at t_phi.

From the renewal relation (their Eq.1),  Qt_r(s) = Qt_0(s+r)/(1 - r Qt_0(s+r)):
    d/dr Qt_r(s)|_{r=0} = Qt_0'(s) + Qt_0(s)^2
and since  L{-t Q_0(t)} = Qt_0'(s)  and  L{(Q_0*Q_0)(t)} = Qt_0(s)^2,

    dF_r(t)/dr|_{r=0}  =  t Q_0(t) - (Q_0 * Q_0)(t)  =:  g(t)

with  Q_0(t) = erf(d / sqrt(4 D t)).  Then t_c solves g(t_c)=0 and

    phi_c = 1 - Q_0(t_c).

Everything below is a plain quadrature of a closed-form function.
"""
import numpy as np
from scipy.special import erf
from scipy.integrate import quad
from scipy.optimize import brentq

def Q0(t, d=1.0, D=1.0):
    return erf(d/np.sqrt(4*D*t)) if t > 0 else 1.0

def conv(t, d=1.0, D=1.0):
    """(Q_0 * Q_0)(t).  Substituting u = t*x removes the endpoint kinks."""
    f = lambda x: Q0(t*x, d, D)*Q0(t*(1-x), d, D)
    v, err = quad(f, 0, 1, limit=400, epsabs=1e-13, epsrel=1e-13)
    return t*v, t*err

def g(t, d=1.0, D=1.0):
    c, _ = conv(t, d, D)
    return t*Q0(t, d, D) - c

if __name__ == "__main__":
    print("g(t) = t Q_0(t) - (Q_0*Q_0)(t)      [resetting helps the quantile iff g>0]")
    for t in [1e-3, 1e-2, 0.05, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 10, 100]:
        print(f"  t={t:8.3g}   g={g(t):+.6e}   Q_0={Q0(t):.6f}")

    # g underflows to 0 at small t (resetting is irrelevant before the target
    # is reachable), so brentq must be bracketed on the REAL crossing, not on
    # the numerically-flat left tail.  Bracket found by scan.
    ts = np.geomspace(1e-2, 1e3, 400); gs = np.array([g(t) for t in ts])
    i = np.where((gs[:-1] < 0) & (gs[1:] > 0))[0]
    assert len(i) == 1, f"expected exactly one sign change, got {len(i)}"
    tc = brentq(g, ts[i[0]], ts[i[0]+1], xtol=1e-14, rtol=8.9e-16)
    phic = 1 - Q0(tc)
    print()
    print(f"  t_c   = {tc:.10f}   (in units d^2/D)")
    print(f"  Q_0   = {Q0(tc):.10f}")
    print(f"  phi_c = {phic:.10f}        paper: 0.412")

    # the number must not depend on d or D -- it is a property of the shape only
    print("\n  scale invariance check (phi_c must be identical):")
    for d, D in [(1,1), (3.7,1), (1,0.21), (2.5,11.0)]:
        sc = d*d/D
        t = brentq(lambda t: g(t,d,D), 0.3*sc, 1.2*sc, xtol=1e-16, rtol=8.9e-16)
        print(f"    d={d:5}, D={D:6}  ->  t_c={t:.6g}   phi_c={1-Q0(t,d,D):.10f}")
