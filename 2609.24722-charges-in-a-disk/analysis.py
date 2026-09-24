"""The decisive test: which basis, fitted to small N, predicts N = 10^5?

Amore & Zarate fit {N^2, N^{3/2}, N, N^{1/2}, 1}.  fitbasis.py argues the
next term after N^{3/2} is N^{4/3}, not N.  Fit my own ladder both ways with
k1 pinned to pi/4, then extrapolate FIFTY TIMES beyond the fitting range and
compare with Lavrov & Nikonov's 31-hour number, which neither fit has seen.
"""
import json
import numpy as np
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


if __name__ == "__main__":
    import sys
    main(*sys.argv[1:])
