#!/usr/bin/env python3
"""exponent.py -- the local exponent of kappa(1-2delta), pushed until it stops
moving, and checked for N-dependence before any number is quoted.

curvature.py established gain = kappa eps^2 with the exponent 2.000 at every
delta once eps is scaled to (1-2delta).  kappa itself diverges as the directed
ring approaches the undirected (Hermitian) point where Theorem 4 does not apply
and the effect is exactly zero.  A single global power-law fit gave -1.40 and a
fit restricted to the tail gave -1.35, which means the fit had not converged and
neither number should be quoted.  This computes the LOCAL exponent between
consecutive points instead, so a drifting slope cannot hide inside an average.
"""
from __future__ import annotations

import numpy as np

from directed_ring import L_directed_ring, best_homogeneous_L, best_perturbation


def kappa_at(N, a, epsrel=(0.1, 0.05, 0.02, 0.01)):
    delta = (1 - a) / 2
    L = L_directed_ring(N, delta)
    bh, Lh = best_homogeneous_L(L)
    g, e = [], []
    for r in epsrel:
        eps = r * a
        base, best, _ = best_perturbation(L, bh, eps)
        d = base - best
        if d > 1e-13:
            g.append(d); e.append(eps)
    if len(g) < 3:
        return np.nan, np.nan, Lh
    p = np.polyfit(np.log(e), np.log(g), 1)
    return float(np.exp(p[1])), float(p[0]), Lh


if __name__ == "__main__":
    AS = [0.02, 0.01, 0.005, 0.002, 0.001, 0.0005, 0.0002, 0.0001]
    for N in (4, 6, 8):
        print(f"\n=== N = {N} ===")
        print(f"{'1-2delta':>10} {'kappa':>13} {'p(eps)':>7} {'|Lmax|':>9} "
              f"{'local exponent':>15}")
        prev = None
        for a in AS:
            k, p, Lh = kappa_at(N, a)
            loc = ""
            if prev and np.isfinite(k):
                loc = f"{np.log(k/prev[1]) / np.log(a/prev[0]):15.4f}"
            print(f"{a:10.5f} {k:13.4f} {p:7.3f} {abs(Lh):9.5f} {loc}")
            if np.isfinite(k):
                prev = (a, k)
    print("\nlocal exponent = d log kappa / d log (1-2delta) between adjacent rows.")
    print("-4/3 = -1.3333.  If all three N agree and the drift is settling there,")
    print("that is a number; if they disagree, the exponent is N-dependent.")
