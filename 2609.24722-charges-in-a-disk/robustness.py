"""Try to kill the out-of-sample result before believing it.

analysis.py says: pin k1 = pi/4 and k2 = -1.5642653, fit only k_{4/3} and k3
on 100 <= N <= 660, and the prediction at N = 10^5 -- 152x beyond the data --
lands 4.3e-7 from the measured value.  Two free parameters.  That is as good
as the published five-parameter fit over a range 7.6x wider, which is exactly
the kind of result I want to be true, so it gets attacked rather than
reported.

Every knob I am unsure about, moved: which points are used, whether k2 is
pinned at theory or at the fitted value, and whether the N^{4/3} term is
there at all.
"""
import json
import numpy as np

K1 = np.pi / 4
K2_TH = -1.5642652795
K2_AZ = -1.5628
LN_N, LN_E = 100_000, 7.80466624157e9


def predict(N, E, powers, k2):
    N = np.asarray(N, float); E = np.asarray(E, float)
    y = E - K1 * N ** 2 - k2 * N ** 1.5
    A = np.stack([N ** p for p in powers], 1)
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    rms = np.sqrt(((y - A @ c) ** 2).mean())
    E5 = K1 * LN_N ** 2 + k2 * LN_N ** 1.5 + sum(ci * LN_N ** p
                                                 for ci, p in zip(c, powers))
    return abs(E5 - LN_E) / LN_E, c, rms


def main(path="ladder.json"):
    rows = json.load(open(path))
    N = np.array([r["N"] for r in rows], float)
    E = np.array([r["E"] for r in rows], float)

    print("=" * 78)
    print("ADVERSARIAL: relative error at N = 1e5, fitted only on the subset shown")
    print("=" * 78)
    print(f"{'subset':<22}{'k2 pinned at':<16}{'basis':<16}{'rel @1e5':>11}{'k_4/3':>10}")
    print("-" * 78)
    subsets = [("all " + str(len(N)), slice(None)),
               ("drop smallest", slice(1, None)),
               ("drop two smallest", slice(2, None)),
               ("drop largest", slice(0, -1)),
               ("drop largest two", slice(0, -2)),
               ("middle only", slice(1, -1))]
    for name, sl in subsets:
        if len(N[sl]) < 3:
            continue
        for k2name, k2 in (("theory", K2_TH), ("Amore-Zarate", K2_AZ)):
            for bname, powers in (("N^{4/3} + N", [4 / 3, 1.0]), ("N only", [1.0])):
                if len(N[sl]) <= len(powers):
                    continue
                rel, c, rms = predict(N[sl], E[sl], powers, k2)
                k43 = c[0] if len(powers) == 2 else float("nan")
                print(f"{name:<22}{k2name:<16}{bname:<16}{rel:11.2e}{k43:10.5f}")
        print("-" * 78)
    print(f"{'published 5-param fit':<54}{2.525e-7:11.2e}")
    print("\n  If the 4e-7 only appears for one subset it is luck. If it appears")
    print("  across subsets with k2 pinned at THEORY and not with k2 pinned at")
    print("  the fitted value, that is evidence the theoretical k2 is the right")
    print("  one and the fitted digits were absorbing the missing term.")


if __name__ == "__main__":
    import sys
    main(*sys.argv[1:])
