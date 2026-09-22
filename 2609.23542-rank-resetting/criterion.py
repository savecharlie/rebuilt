"""
One criterion containing both limits of the paper.

<T_(k)>(r) = int_0^inf S_k(Q_r(t)) dt,  S_k(Q) = sum_{m<k} C(N,m)(1-Q)^m Q^{N-m}.
Differentiating at r=0 and using  dQ_r/dr|_0 = -g(t),  g = t Q_0 - (Q_0*Q_0):

    d<T_(k)>/dr |_{r=0}  =  - int_0^inf S_k'(Q_0) g(t) dt

and the standard binomial-CDF identity  S_k'(Q) = N C(N-1,k-1) (1-Q)^{k-1} Q^{N-k}
(positive constant) leaves the sign entirely to

    J(N,k) = int_0^inf (1 - Q_0(t))^{k-1} Q_0(t)^{N-k} g(t) dt .

    J > 0  <=>  a finite optimal reset rate exists (resetting helps at r=0+).

Tail: integrand ~ t^{(1-N+k)/2}, so J = +inf for k >= N-3 -- resetting ALWAYS
helps those ranks.  The criterion only bites for k <= N-4.

Limits:
  k=1              -> critical population N_c for the fastest arrival
  N->inf, k/N=phi  -> weight concentrates on Q_0 = 1-phi, giving sign(g(t_phi)),
                      i.e. exactly the phi_c threshold.
"""
import numpy as np
from scipy.special import erf, erfc, comb
from scipy.integrate import quad
from scipy.optimize import brentq
from closedform import g_closed

def Q0(t): return erf(1/(2*np.sqrt(t)))

def J(N, k, lo=-12.0, hi=None):
    """int_0^inf (1-Q0)^{k-1} Q0^{N-k} g dt, over t = e^u."""
    if hi is None:
        hi = 6.0 + 2.2*max(N-k, 1)           # far into the algebraic tail
    f = lambda u: ((1-Q0(np.exp(u)))**(k-1) * Q0(np.exp(u))**(N-k)
                   * g_closed(np.exp(u)) * np.exp(u))
    v, e = quad(f, lo, hi, limit=600, epsabs=1e-14, epsrel=1e-12)
    return v, e

def Sk_prime_check():
    """numerically verify S_k'(Q) = N C(N-1,k-1)(1-Q)^{k-1} Q^{N-k}"""
    rng = np.random.default_rng(1); worst = 0.0
    for N in range(1, 12):
        for k in range(1, N+1):
            for Q in rng.uniform(0.05, 0.95, 5):
                h = 1e-6
                S = lambda q: sum(comb(N,m,exact=True)*(1-q)**m*q**(N-m) for m in range(k))
                num = (S(Q+h)-S(Q-h))/(2*h)
                ana = N*comb(N-1,k-1,exact=True)*(1-Q)**(k-1)*Q**(N-k)
                worst = max(worst, abs(num-ana)/max(abs(ana),1e-12))
    return worst

if __name__ == "__main__":
    print(f"binomial-CDF derivative identity, worst relative error: {Sk_prime_check():.2e}")
    print()
    print("J(N,1) -- the fastest arrival.  Sign change = critical population N_c.")
    print(f"{'N':>5} {'J(N,1)':>16}")
    for N in [5, 6, 7, 7.2, 7.3, 7.33, 7.4, 8, 10]:
        v, e = J(N, 1); print(f"{N:5g} {v:+16.9f}")
    Nc = brentq(lambda N: J(N,1)[0], 5.0, 12.0, xtol=1e-10)
    print(f"\n  N_c (fastest arrival, k=1) = {Nc:.8f}")
    print()
    print("N=6, all ranks.  Compare with the direct optimisation: r*_k > 0 for every k.")
    for k in range(1, 7):
        if k >= 6-3:
            print(f"   k={k}:  J = +inf  (tail exponent {(1-6+k)/2:+.1f}; resetting always helps)")
        else:
            v,e = J(6,k); print(f"   k={k}:  J = {v:+.8f}")
