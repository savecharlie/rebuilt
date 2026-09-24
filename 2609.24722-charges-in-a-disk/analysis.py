"""The decisive test: which basis, fitted to small N, predicts N = 10^5?

Amore & Zarate fit {N^2, N^{3/2}, N, N^{1/2}, 1}.  fitbasis.py argues the
next term after N^{3/2} is N^{4/3}, not N.  Fit my own ladder both ways with
k1 pinned to pi/4, then extrapolate FIFTY TIMES beyond the fitting range and
compare with Lavrov & Nikonov's 31-hour number, which neither fit has seen.
"""
import json, sys
import numpy as np
sys.path.insert(0, "/home/ivy/iris-the-maker/tools")
import dull
from fitbasis import fit

K1 = np.pi / 4
K2_THEORY = -1.5642652795
LN_N, LN_E = 100_000, 7.80466624157e9
AZ = {1.5: -1.5628, 1.0: 1.0302, 0.5: -0.9899, 0.0: 4.9255}

BASES = {
    "Amore-Zarate {3/2, 1, 1/2, 0}":        [2.0, 1.5, 1.0, 0.5, 0.0],
    "with N^{4/3}  {3/2, 4/3, 1, 1/2}":     [2.0, 1.5, 4 / 3, 1.0, 0.5],
    "two terms only {3/2}":                 [2.0, 1.5],
}


def main(path="ladder.json"):
    rows = json.load(open(path))
    N = np.array([r["N"] for r in rows], float)
    E = np.array([r["E"] for r in rows], float)
    print(f"  ladder: N = {[int(n) for n in N]}")
    print(f"  (local minima sit ABOVE the true ones, which biases k2 toward")
    print(f"   SMALLER magnitude — i.e. against the claim being tested here.)\n")

    print("=" * 78)
    print(f"{'basis':<36}{'k2':>12}{'rms resid':>12}{'pred N=1e5':>16}{'rel':>10}")
    print("=" * 78)
    for name, powers in BASES.items():
        c, res = fit(N, E, powers, pin={2.0: K1})
        pred = sum(c[p] * LN_N ** p for p in powers)
        print(f"{name:<36}{c[1.5]:12.6f}{np.sqrt((res**2).mean()):12.3e}"
              f"{pred:16.7e}{abs(pred-LN_E)/LN_E:10.2e}")
        # same basis, k2 PINNED at theory
        c2, res2 = fit(N, E, powers, pin={2.0: K1, 1.5: K2_THEORY})
        pred2 = sum(c2[p] * LN_N ** p for p in powers)
        print(f"{'   ... same, k2 pinned at theory':<36}{K2_THEORY:12.6f}"
              f"{np.sqrt((res2**2).mean()):12.3e}{pred2:16.7e}"
              f"{abs(pred2-LN_E)/LN_E:10.2e}")
    print("=" * 78)
    print(f"{'Lavrov & Nikonov, measured':<36}{'':12}{'':12}{LN_E:16.7e}")
    print(f"{'Amore-Zarate published fit':<36}{AZ[1.5]:12.6f}{'':12}"
          f"{K1*LN_N**2 + sum(AZ[p]*LN_N**p for p in AZ):16.7e}"
          f"{abs(K1*LN_N**2 + sum(AZ[p]*LN_N**p for p in AZ)-LN_E)/LN_E:10.2e}")

    print("\n" + "=" * 78)
    print("SLIDING WINDOW: does the fitted k2 move toward theory as N grows?")
    print("=" * 78)
    powers = BASES["Amore-Zarate {3/2, 1, 1/2, 0}"]
    for lo in range(0, max(1, len(N) - 4)):
        sel = slice(lo, None)
        c, _ = fit(N[sel], E[sel], powers, pin={2.0: K1})
        print(f"  N >= {int(N[lo]):5d}   ({len(N[sel])} points)   "
              f"k2 = {c[1.5]:.6f}   gap to theory {c[1.5]-K2_THEORY:+.6f}")
    resolving_power(N, E)
    deficit_fit(N, E)


def resolving_power(N, E):
    """Is -1.5628 even distinguishable from -1.5642653 by a fit like theirs?

    Fit my own ladder in the Amore-Zarate basis with k1 pinned, and carry the
    least-squares covariance through.  If the standard error on k2 is bigger
    than the 0.0015 gap, then the fitted value was never determined to the
    precision its four digits imply, and there is nothing to explain.
    """
    print("\n" + "=" * 78)
    print("CAN A FIT LIKE THEIRS EVEN SEE THE DIFFERENCE?")
    print("=" * 78)
    N = np.asarray(N, float); E = np.asarray(E, float)
    powers = [1.5, 1.0, 0.5, 0.0]
    A = np.stack([N ** p for p in powers], 1)
    y = E - K1 * N ** 2
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = y - A @ c
    dof = len(N) - len(powers)
    s2 = float(resid @ resid) / max(dof, 1)
    cov = s2 * np.linalg.inv(A.T @ A)
    se = np.sqrt(np.diag(cov))
    # condition number after scaling each column to unit norm
    An = A / np.linalg.norm(A, axis=0)
    print(f"  design-matrix condition number (columns normalised): "
          f"{np.linalg.cond(An):.3e}")
    for p, ci, si in zip(powers, c, se):
        print(f"     k_{{{p:.3g}}} = {ci:+.6f}  +- {si:.6f}")
    k2, sk2 = c[0], se[0]
    print(f"\n  fitted k2      {k2:+.6f} +- {sk2:.6f}")
    print(f"  theory         {K2_THEORY:+.6f}   ({abs(k2-K2_THEORY)/sk2:.2f} sigma away)")
    print(f"  Amore-Zarate   {AZ[1.5]:+.6f}   ({abs(AZ[1.5]-K2_THEORY)/sk2:.2f} sigma"
          f" from theory on this scale)")
    return k2, sk2


def deficit_fit(N, E):
    """The cleanest statement available: pin BOTH theorems, fit what is left.

        k2_eff(N) = k2 - [ k_{4/3} N^{-1/6} + k3 N^{-1/2} + k4 N^{-1} + ... ]

    so the residual after removing pi/4 and the Madelung term is a pure
    boundary series with only two or three numbers in it.  Fit it on my own
    ladder alone and then ask what it says at N = 10^5, which is fifty times
    past the end of the data and which the fit has never seen.
    """
    print("\n" + "=" * 78)
    print("PIN BOTH THEOREMS, FIT ONLY THE BOUNDARY SERIES")
    print("=" * 78)
    N = np.asarray(N, float); E = np.asarray(E, float)
    keff = -(K1 * N ** 2 - E) / N ** 1.5
    deficit = keff - K2_THEORY          # positive, shrinking:
    #   k2_eff(N) = k2 + k_{4/3} N^{-1/6} + k3 N^{-1/2} + k4 N^{-1} + ...
    for name, exps in (("N^{4/3} + N", (-1/6, -0.5)),
                       ("N only (Amore-Zarate basis)", (-0.5,)),
                       ("N^{4/3} + N + N^{1/2}", (-1/6, -0.5, -1.0))):
        A = np.stack([N ** e for e in exps], 1)
        c, *_ = np.linalg.lstsq(A, deficit, rcond=None)
        pred_def = sum(ci * N ** e for ci, e in zip(c, exps))
        rms = np.sqrt(((deficit - pred_def) ** 2).mean())
        # out of sample
        d5 = sum(ci * LN_N ** e for ci, e in zip(c, exps))
        E5 = K1 * LN_N ** 2 + (K2_THEORY + d5) * LN_N ** 1.5
        coef = "  ".join(f"k_{{{1.5+e:.3g}}}={ci:+.5f}" for ci, e in zip(c, exps))
        print(f"  {name:<28} rms {rms:.2e}   {coef}")
        print(f"  {'':<28} -> E(1e5) = {E5:.7e}   rel to measured "
              f"{abs(E5-LN_E)/LN_E:.2e}")
    print(f"  {'measured (Lavrov-Nikonov)':<28} -> E(1e5) = {LN_E:.7e}")
    print(f"\n  running coefficient, ladder then their datum:")
    for n, k in zip(N, keff):
        print(f"     N={int(n):6d}  k2_eff = {k:.6f}   gap to theory "
              f"{k-K2_THEORY:+.6f}")
    k5 = -(K1 * LN_N ** 2 - LN_E) / LN_N ** 1.5
    print(f"     N={LN_N:6d}  k2_eff = {k5:.6f}   gap to theory "
          f"{k5-K2_THEORY:+.6f}   <- 31 CPU-hours")
    print(f"     theory       {K2_THEORY:.6f}")

    # MEASURE the exponent of the deficit rather than writing one in prose.
    # N^{4/3} in the energy  ->  deficit ~ N^{-1/6}
    # N      in the energy  ->  deficit ~ N^{-1/2}
    allN = np.append(N, LN_N)
    alld = np.append(deficit, k5 - K2_THEORY)
    print("\n  deficit exponent, measured (dull.scaling):")
    with dull.watch() as w:
        p_all = dull.scaling(allN, alld, "deficit, 100..1e5")
        p_lad = dull.scaling(N, deficit, "deficit, ladder alone")
        print(f"     100 .. 1e5  (ladder + their datum): {p_all:+.4f}")
        print(f"     100 .. {int(N[-1]):<4d}  (ladder alone)      : {p_lad:+.4f}")
        print(f"     pure N^(4/3) term would give -1/6 = {-1/6:+.4f}")
        print(f"     pure N       term would give -1/2 = {-0.5:+.4f}")
        w.report()


if __name__ == "__main__":
    import sys
    main(*sys.argv[1:])
