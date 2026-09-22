"""
The bridge: one integral, both published limits, and the finite-N correction.

    J(N,k) = int_0^inf (1-Q_0)^{k-1} Q_0^{N-k} g(t) dt

  k = 1          ->  sign change at N_c = 7.32647733   (Biroli/Majumdar/Schehr 2023: 7.3264773...)
  N -> inf, k/N  ->  sign change at phi_c = 0.41231018  (Belan 2020: 0.4123)

NUMERICS.  The weight (1-Q)^{k-1} Q^{N-k} underflows float64 above N ~ 1000
(Q^755 ~ e^-400).  A first version of this file did NOT do that and silently
returned phi_c = 0.4995 for N = 1280, 2560 -- the integrand had gone to zero
everywhere and brentq found a spurious root at k = N/2.  The weight is now
carried in logs with its peak factored out, which costs nothing and cannot
underflow.  Sign of J is unaffected by the positive constant removed.
"""
import numpy as np
from scipy.special import erf
from scipy.integrate import quad
from scipy.optimize import brentq
from closedform import g_closed

def Q0(t): return erf(1/(2*np.sqrt(t)))
def t_of_Q(q): return brentq(lambda lt: Q0(np.exp(lt)) - q, -12, 20, xtol=1e-13)

def J(N, k):
    Qpk = min(max((N-k)/(N-1.0), 1e-12), 1-1e-12)
    lt0 = t_of_Q(Qpk)
    logpk = (k-1)*np.log1p(-Qpk) + (N-k)*np.log(Qpk)   # peak of the log-weight
    w = max(6.0/np.sqrt(N), 0.05)
    def f(u):
        t = np.exp(u); q = Q0(t)
        lw = (k-1)*np.log1p(-q) + (N-k)*np.log(q) - logpk
        return np.exp(np.clip(lw, -700, 700)) * g_closed(t) * t
    edges = sorted({max(-14.0, lt0-60*w), lt0-14*w, lt0-4*w, lt0-w,
                    lt0+w, lt0+4*w, lt0+14*w, lt0+60*w, lt0+200*w})
    return sum(quad(f, a, b, limit=400, epsabs=1e-300, epsrel=1e-11)[0]
               for a, b in zip(edges[:-1], edges[1:]))

if __name__ == "__main__":
    PHI = 0.412310175459
    print("phi_c(N) = k_c/N, where J(N,k_c)=0.        asymptote: 0.412310175")
    print(f"{'N':>7} {'k_c':>13} {'phi_c(N)':>14} {'deviation':>13} {'dev*N':>9}")
    for N in [12, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240]:
        kc = brentq(lambda k: J(N, k), 1.0+1e-9, N-3.5, xtol=1e-10)
        phi = kc/N; dev = phi - PHI
        print(f"{N:7d} {kc:13.5f} {phi:14.9f} {dev:+13.3e} {dev*N:+9.4f}")
    print()
    print("  dev * N is constant -> phi_c(N) = phi_c - c/N, a 1/N approach, not 1/sqrt(N).")
    print("  (An earlier run of this file printed dev*sqrt(N) and I wrote the wrong law")
    print("   under it; the numbers had said 1/N all along -- each doubling halved dev.)")
