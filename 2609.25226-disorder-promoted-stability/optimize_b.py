#!/usr/bin/env python3
"""optimize_b.py -- does heterogeneous damping beat the best homogeneous damping?

Runs only against instruments `stability.py --selftest` has already checked
against the closed form.  The homogeneous benchmark is NOT scanned here; it is
the exact -sqrt(mu_2), so the comparison has no numerical optimum on the side it
is trying to beat.  Any gain reported below therefore cannot be a search
artefact in the baseline.

Lambda_max is kinked (eigenvalues collide and split into complex pairs at exactly
the damping values that matter), so Nelder-Mead from many starts, not a gradient.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize

from stability import J_second, J_leaky, lam_max, laplacian


def optimise(model, A, starts=40, seed=0, bmax=50.0, budget=None,
             seed_at=None, **kw):
    """min_b Lambda_max.  `budget` fixes sum(b) (the 1D case needs it: there,
    more damping is always better and the unconstrained optimum runs to the
    bound, where every configuration is trivially at the corner).

    SEED AT THE HOMOGENEOUS OPTIMUM.  The first version of this did not, and on
    `complete N=6` it reported the heterogeneous optimum as 23% WORSE than the
    homogeneous one -- which is impossible, since the homogeneous point is in
    the feasible set.  That was Nelder-Mead failing to reach a kink from a
    random start, and I nearly wrote it down as evidence that heterogeneity
    does not help on undirected networks.  With `seed_at` in the start list the
    optimiser can only tie or win, so a zero is a real local statement and a
    negative number is impossible by construction."""
    N = A.shape[0]
    rng = np.random.default_rng(seed)

    def unpack(z):
        b = np.abs(z)
        if budget is not None:
            s = b.sum()
            b = b * (budget / s) if s > 1e-12 else np.full(N, budget / N)
        return np.clip(b, 0.0, bmax)

    def obj(z):
        return lam_max(model(unpack(z), A, **kw))

    base = np.full(N, seed_at) if seed_at is not None else np.full(N, 1.0)
    best = (obj(base), unpack(base))
    for t in range(starts):
        z0 = base * rng.uniform(0.3, 2.0, N) if t else base.copy()
        r = minimize(obj, z0, method="Nelder-Mead",
                     options={"maxiter": 20000, "xatol": 1e-10, "fatol": 1e-12})
        if r.fun < best[0]:
            best = (float(r.fun), unpack(r.x))
    return best


def spread(b):
    return float(np.std(b) / max(1e-12, np.mean(b)))


def ring(N, k=1.0):
    A = np.zeros((N, N))
    for i in range(N):
        A[i, (i + 1) % N] = A[(i + 1) % N, i] = k
    return A


def complete(N, k=1.0):
    A = k * (np.ones((N, N)) - np.eye(N))
    return A


if __name__ == "__main__":
    print("=== 2D: second-order Kuramoto, J = [[0,I],[-L,-B]] ===")
    print("homogeneous benchmark is the EXACT -sqrt(mu_2), not a scan.\n")
    print(f"{'network':16} {'mu_2':>8} {'hom Lmax':>10} {'het Lmax':>10} "
          f"{'gain':>8}  {'spread(b*)':>10}")
    for name, A in [("ring N=4", ring(4)), ("ring N=6", ring(6)),
                    ("ring N=8", ring(8)), ("complete N=4", complete(4)),
                    ("complete N=6", complete(6)), ("path N=5", None)]:
        if A is None:
            A = np.zeros((5, 5))
            for i in range(4):
                A[i, i + 1] = A[i + 1, i] = 1.0
        mu = np.linalg.eigvalsh(laplacian(A))
        hom = -np.sqrt(mu[1])
        het, bstar = optimise(J_second, A, starts=30, seed=1,
                              seed_at=2 * np.sqrt(mu[1]))
        print(f"{name:16} {mu[1]:8.4f} {hom:10.5f} {het:10.5f} "
              f"{(hom-het)/abs(hom)*100:7.2f}%  {spread(bstar):10.4f}")
        print(f"{'':16} b* = " + " ".join(f"{v:.4f}" for v in bstar))

    print("\n=== 1D control: leaky Kuramoto, J = -(B+L), Hermitian ===")
    print("budget sum(b) fixed, so the comparison is not decided by the bound.\n")
    print(f"{'network':16} {'budget':>7} {'hom Lmax':>10} {'het Lmax':>10} "
          f"{'gain':>8}  {'spread(b*)':>10}")
    for name, A in [("ring N=6", ring(6)), ("complete N=5", complete(5)),
                    ("path N=5", None)]:
        if A is None:
            A = np.zeros((5, 5))
            for i in range(4):
                A[i, i + 1] = A[i + 1, i] = 1.0
        N = A.shape[0]
        for budget in (N * 1.0,):
            hom = lam_max(J_leaky(np.full(N, budget / N), A), drop_null=False)
            het, bstar = optimise(J_leaky, A, starts=30, seed=2, budget=budget,
                                  seed_at=budget / N)
            print(f"{name:16} {budget:7.2f} {hom:10.5f} {het:10.5f} "
                  f"{(hom-het)/abs(hom)*100:7.2f}%  {spread(bstar):10.4f}")
            print(f"{'':16} b* = " + " ".join(f"{v:.4f}" for v in bstar))
