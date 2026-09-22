"""
Derive the 1/N coefficient instead of fitting it.

Change variables Q = Q_0(t).  Since Q_0(t) = erf(1/(2 sqrt t)),
    T(Q) = 1/(4 y^2),   |T'(Q)| = (sqrt(pi)/4) e^{y^2} / y^3,   y = erfinv(Q).
Then
    J(N,k) = int_0^1 (1-Q)^{k-1} Q^{N-k} psi(Q) dQ,   psi(Q) = g(T(Q)) |T'(Q)|,
so J is proportional to E[psi(Q)] under Beta(a=N-k+1, b=k):
    mean   mu    = (N-k+1)/(N+1) = (1-phi) + phi/(N+1)
    var    s^2   = ab/((N+1)^2 (N+2)) ~ phi(1-phi)/N
Laplace expansion:  E[psi] = psi(mu) + (1/2) psi''(mu) s^2 + O(1/N^2).
J = 0 with psi(Q*) = 0 at Q* = 1 - phi_c gives

    phi_c(N) = phi_c - c/N,     c = phi_c + (1/2)(psi''/psi')(Q*) phi_c (1-phi_c)
"""
import numpy as np
from scipy.special import erf, erfinv
from closedform import g_closed

PHI = 0.412310175459
QSTAR = 1 - PHI

def psi(Q):
    y = erfinv(np.asarray(Q, dtype=float))
    return g_closed(1.0/(4*y*y)) * (np.sqrt(np.pi)/4) * np.exp(y*y)/y**3

def deriv(f, x, n, h):
    if n == 1: return (f(x+h)-f(x-h))/(2*h)
    if n == 2: return (f(x+h)-2*f(x)+f(x-h))/(h*h)

print("psi(Q*) must vanish (it is g at its own root):", f"{psi(QSTAR):.3e}")
print()
print("%10s %14s %14s %12s %12s" % ("h","psi1","psi2","ratio","c"))
for h in [2e-2, 1e-2, 5e-3, 2.5e-3, 1.25e-3]:
    d1 = deriv(psi, QSTAR, 1, h); d2 = deriv(psi, QSTAR, 2, h)
    c = PHI + 0.5*(d2/d1)*PHI*(1-PHI)
    print(f"{h:10.2e} {d1:14.6f} {d2:14.6f} {d2/d1:12.5f} {c:12.6f}")
print()
print("  measured from the N-sweep:  dev*N -> -2.0763,  i.e.  c = +2.0763")
