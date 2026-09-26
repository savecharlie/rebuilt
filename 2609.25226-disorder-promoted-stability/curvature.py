#!/usr/bin/env python3
"""curvature.py -- HOW BIG is disorder-promoted stability, as the network
approaches the Hermitian point where it provably vanishes?

Theorem 4 proves the effect EXISTS for every directed first-neighbour ring
(delta != 1/2) and Proposition 5 says the homogeneous optimum is a kink on
undirected ones.  Neither says what happens in between.  The directed_ring.py
run found the gain is exactly quadratic in the perturbation, gain = kappa eps^2,
with kappa jumping around by two and a half orders of magnitude across delta.
This measures kappa(delta) properly: at each delta, fit gain against eps over a
decade and check the exponent really is 2 before quoting a curvature.

An exponent that comes back at 2 is the negative Hessian eigenvalue of Theorem 4.
An exponent near 1 would be a kink, and kappa would be meaningless.  The fit
reports the exponent every time so a bad one cannot hide inside a number.
"""
from __future__ import annotations

import numpy as np

from directed_ring import L_directed_ring, J_of, best_homogeneous_L, best_perturbation


def kappa(L, bh, eps_list):
    g, e = [], []
    for eps in eps_list:
        base, best, pair = best_perturbation(L, bh, eps)
        d = base - best
        if d > 1e-14:
            g.append(d); e.append(eps)
    if len(g) < 3:
        return np.nan, np.nan, 0
    p = np.polyfit(np.log(e), np.log(g), 1)
    return float(np.exp(p[1])), float(p[0]), len(g)


if __name__ == "__main__":
    N = 6
    # EPS MUST SCALE WITH (1-2delta).  A fixed list looked fine down to
    # 1-2delta = 0.02 and then returned a fitted exponent of 1.67 at 0.002 --
    # not a kink appearing, just eps leaving the quadratic regime, because the
    # regime shrinks with the distance from the Hermitian point.  A curvature
    # quoted off that row would have been a measurement of my own step size.
    EPSREL = [0.2, 0.1, 0.05, 0.02, 0.01]
    print(f"directed ring N={N}.  gain = kappa * eps^p, eps fitted over "
          f"{EPSREL} x (1-2delta)\n")
    print(f"{'delta':>7} {'1-2d':>8} {'|Lmax|':>9} {'kappa':>12} {'p':>7} "
          f"{'kappa(1-2d)':>12} {'kappa(1-2d)^2':>14}")
    rows = []
    for delta in [0.0, 0.05, 0.1, 0.2, 0.3, 0.35, 0.4, 0.44, 0.46, 0.47, 0.48,
                  0.49, 0.495, 0.498, 0.499, 0.4995, 0.4998]:
        L = L_directed_ring(N, delta)
        bh, Lh = best_homogeneous_L(L)
        k, p, n = kappa(L, bh, [e * (1 - 2 * delta) for e in EPSREL])
        a = 1 - 2 * delta
        rows.append((delta, a, Lh, k, p))
        print(f"{delta:7.3f} {a:8.4f} {abs(Lh):9.5f} {k:12.5f} {p:7.3f} "
              f"{k*a:12.5f} {k*a*a:14.5f}")

    # the Hermitian point itself
    L = L_directed_ring(N, 0.5)
    bh, Lh = best_homogeneous_L(L)
    base, best, pair = best_perturbation(L, bh, 0.001)
    print(f"{0.5:7.3f} {0.0:8.4f} {abs(Lh):9.5f}   gain at eps=1e-3: "
          f"{base-best:.3e}  ({'none' if pair is None else pair})")

    good = [(a, k) for (_, a, _, k, p) in rows if abs(p - 2) < 0.02 and a > 0]
    print(f"\n{len(good)} of {len(rows)} rows have a fitted exponent within 0.02 of 2.")
    A = np.array([g[0] for g in good]); K = np.array([g[1] for g in good])
    print(f"{'cut':>8} {'n':>3} {'exponent':>10} {'prefactor':>11} {'max resid':>10}")
    for cut in (1.01, 0.3, 0.1, 0.03, 0.01):
        m = A <= cut
        if m.sum() < 3:
            continue
        f = np.polyfit(np.log(A[m]), np.log(K[m]), 1)
        pr = np.exp(np.polyval(f, np.log(A[m])))
        print(f"{cut:8.3f} {int(m.sum()):3} {f[0]:10.4f} {np.exp(f[1]):11.5f} "
              f"{np.max(np.abs(pr / K[m] - 1)):9.2%}")
    print("If the exponent is still moving at the smallest cut, it has not "
          "converged and no single power should be quoted.")
