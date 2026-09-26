#!/usr/bin/env python3
"""stability.py -- the Jacobians of Montanari, Zanin & Motter (arXiv:2609.25226),
and Lambda_max, with a closed-form check standing in front of every new number.

The paper's claim, rebuilt from Table S1 and nothing else:

  1D nodal dynamics (leaky Kuramoto, consensus):  J = -(B + L),  HERMITIAN for
    undirected A.  Lambda_max(b) is then a CONVEX function of b, so on any convex
    symmetric feasible set some optimum is homogeneous.  Symmetry survives.

  2D nodal dynamics (second-order Kuramoto, spring-mass):  J = [[0, I], [-L, -B]],
    NON-HERMITIAN.  Lambda_max(b) is not convex and the optimum can break the
    symmetry of the network.

THE KNOWN ANSWER THIS IS VALIDATED AGAINST, derived before any optimiser ran.
With homogeneous damping B = b*I the second-order block matrix decouples in the
eigenbasis of L: each Laplacian eigenvalue mu contributes the roots of

    s^2 + b s + mu = 0      ->      s = (-b +- sqrt(b^2 - 4 mu)) / 2

so a mode is underdamped (Re = -b/2) while b < 2 sqrt(mu) and overdamped after.
Raising b helps every underdamped mode and hurts every overdamped one, and the
FIRST mode to overdamp is the smallest positive eigenvalue, the algebraic
connectivity mu_2.  Therefore the best homogeneous damping is exactly

    b* = 2 sqrt(mu_2)        Lambda_max = -sqrt(mu_2)

That is a hard target with no fitted parameter in it, and `--selftest` checks the
numerical Lambda_max against it on random graphs before anything else is believed.

    python3 stability.py --selftest
"""
from __future__ import annotations

import argparse

import numpy as np


# --------------------------------------------------------------------------- #
# Jacobians, verbatim from Table S1.

def laplacian(A):
    """L = diag(row sums) - A.  Out-degree convention, as in the paper."""
    return np.diag(A.sum(axis=1)) - A


def J_leaky(b, A):
    """Leaky Kuramoto / consensus, 1D nodes.  J = -(B + L).  Hermitian iff A is."""
    return -(np.diag(b) + laplacian(A))


def J_second(b, A):
    """Second-order Kuramoto (power grids), 2D nodes.  J = [[0, I], [-L, -B]]."""
    n = len(b)
    L = laplacian(A)
    return np.block([[np.zeros((n, n)), np.eye(n)],
                     [-L, -np.diag(b)]])


def J_spring(b, A, k=1.0):
    """Spring-mass.  J = [[0, I], [-(L + kI), -B]].  k lifts the null mode."""
    n = len(b)
    L = laplacian(A)
    return np.block([[np.zeros((n, n)), np.eye(n)],
                     [-(L + k * np.eye(n)), -np.diag(b)]])


# --------------------------------------------------------------------------- #

_NULL_TOL = 1e-9


def lam_max(J, drop_null=True):
    """Largest transverse Lyapunov exponent.

    The paper excludes eigenvalues that are identically zero for every b (the set
    Z): in the second-order Kuramoto model the global phase shift is one such
    direction and it never decays, so leaving it in would report Lambda_max = 0
    for every configuration and hide the entire effect.  Exactly ONE is dropped,
    and only if it really is at the origin -- dropping by count would silently
    eat a genuine instability that happened to sit near zero.
    """
    ev = np.linalg.eigvals(J)
    re = np.sort(ev.real)[::-1]
    if drop_null:
        i = int(np.argmin(np.abs(ev)))
        if abs(ev[i]) < _NULL_TOL:
            ev = np.delete(ev, i)
            re = np.sort(ev.real)[::-1]
    return float(re[0])


def best_homogeneous(model, A, lo=1e-4, hi=50.0, n=20001, **kw):
    """Scan b over a homogeneous ray.  Returns (b, Lambda_max).

    A scan and not a solver: Lambda_max(b) is exactly the kinked, non-smooth
    object the paper is about, and a gradient method on it reports whichever
    kink it fell into.
    """
    N = A.shape[0]
    bs = np.linspace(lo, hi, n)
    vals = [lam_max(model(np.full(N, b), A, **kw)) for b in bs]
    i = int(np.argmin(vals))
    return float(bs[i]), float(vals[i])


# --------------------------------------------------------------------------- #

def selftest(seed=0):
    rng = np.random.default_rng(seed)
    ok = True

    def chk(name, cond, extra=""):
        nonlocal ok
        print(("  ok   " if cond else "  FAIL ") + name + (("  " + extra) if extra else ""))
        ok = ok and bool(cond)

    # -- 1. the closed form, mode by mode, on random connected graphs ---------
    worst = 0.0
    for t in range(12):
        N = int(rng.integers(4, 11))
        M = rng.random((N, N)) < 0.5
        A = np.triu(M, 1).astype(float) * rng.uniform(0.3, 2.0, (N, N))
        A = A + A.T
        mu = np.linalg.eigvalsh(laplacian(A))
        if mu[1] < 1e-6:               # disconnected: the closed form does not apply
            continue
        b = float(rng.uniform(0.2, 6.0))
        num = np.sort(np.linalg.eigvals(J_second(np.full(N, b), A)).real)[::-1]
        pred = []
        for m in mu:
            d = complex(b * b - 4 * m) ** 0.5
            pred += [((-b + d) / 2).real, ((-b - d) / 2).real]
        pred = np.sort(np.array(pred))[::-1]
        worst = max(worst, float(np.max(np.abs(num - pred))))
    chk("homogeneous 2nd-order spectrum == roots of s^2+bs+mu", worst < 1e-8,
        f"max |num-pred| = {worst:.2e}")

    # -- 2. the optimum of the homogeneous ray is at 2 sqrt(mu_2) -------------
    bad = 0.0
    badL = 0.0
    for t in range(6):
        N = int(rng.integers(4, 9))
        A = np.triu(rng.random((N, N)) < 0.7, 1).astype(float)
        A = A + A.T
        mu = np.linalg.eigvalsh(laplacian(A))
        if mu[1] < 1e-6:
            continue
        b, L = best_homogeneous(J_second, A, lo=1e-3, hi=12.0, n=24001)
        bad = max(bad, abs(b - 2 * np.sqrt(mu[1])) / (2 * np.sqrt(mu[1])))
        badL = max(badL, abs(L + np.sqrt(mu[1])))
    chk("best homogeneous damping == 2 sqrt(mu_2)", bad < 2e-3, f"rel err {bad:.1e}")
    chk("and it gives Lambda_max == -sqrt(mu_2)", badL < 2e-3, f"abs err {badL:.1e}")

    # -- 3. Hermiticity, which is the whole mechanism -------------------------
    N = 6
    A = np.triu(rng.random((N, N)), 1); A = A + A.T
    b = rng.uniform(0.1, 3, N)
    Jl = J_leaky(b, A)
    chk("1D Jacobian is Hermitian for undirected A", np.allclose(Jl, Jl.T))
    chk("1D spectrum is real", np.max(np.abs(np.linalg.eigvals(Jl).imag)) < 1e-12)
    Js = J_second(b, A)
    chk("2D Jacobian is NOT Hermitian", not np.allclose(Js, Js.T))

    # -- 4. the null direction is real, and dropped exactly once --------------
    J = J_second(np.full(5, 1.3), np.triu(np.ones((5, 5)), 1) + np.triu(np.ones((5, 5)), 1).T)
    allev = np.linalg.eigvals(J)
    chk("2nd-order Kuramoto has exactly one zero eigenvalue",
        int(np.sum(np.abs(allev) < 1e-9)) == 1)
    chk("lam_max drops it", lam_max(J) < -1e-6, f"{lam_max(J):.4f}")
    chk("lam_max WITHOUT dropping it is 0", abs(lam_max(J, drop_null=False)) < 1e-9)

    # -- 5. spring-mass has no null mode, so nothing may be dropped -----------
    A = np.array([[0., 1.], [1., 0.]])
    Jk = J_spring(np.array([1.0, 1.0]), A, k=1.0)
    chk("spring-mass (k>0) has no zero eigenvalue",
        np.min(np.abs(np.linalg.eigvals(Jk))) > 1e-6)
    chk("lam_max leaves it alone",
        abs(lam_max(Jk) - lam_max(Jk, drop_null=False)) < 1e-12)

    print("SELFTEST", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    raise SystemExit(selftest() if a.selftest else selftest())
