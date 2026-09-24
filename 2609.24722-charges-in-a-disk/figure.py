"""The whole claim in one picture.

Plot the RUNNING second coefficient,

        k2_eff(N) = -[ (pi/4) N^2 - E(N) ] / N^{3/2},

which is what k2 would be if the expansion stopped at two terms.  It has to
approach the true k2 from above (in magnitude, from below), and the question
is what it approaches.
"""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

K1 = np.pi / 4
K2_THEORY = -1.5642652795            # -C_M * int rho^{3/2}, see k2.py
AZ = dict(k2=-1.5628, k3=1.0302, k4=-0.9899, k5=4.9255)
LN = (100_000, 7.80466624157e9)


def keff(N, E):
    return -(K1 * np.asarray(N, float) ** 2 - np.asarray(E, float)) / np.asarray(N, float) ** 1.5


def main(path="ladder.json", out="running_k2.png"):
    rows = json.load(open(path))
    N = np.array([r["N"] for r in rows], float)
    E = np.array([r["E"] for r in rows], float)

    fig, ax = plt.subplots(figsize=(8.2, 5.4), dpi=170)
    Ng = np.geomspace(60, 3e5, 500)
    az = (AZ["k2"] * Ng ** 1.5 + AZ["k3"] * Ng + AZ["k4"] * np.sqrt(Ng) + AZ["k5"])
    ax.plot(Ng, az / Ng ** 1.5, color="#c0392b", lw=1.6,
            label="Amore–Zárate 5-term fit (100 ≤ N ≤ 5000)")

    ax.axhline(K2_THEORY, color="#1a6b54", lw=2.2, ls="-",
               label=r"$k_2=-C_M\!\int\!\rho^{3/2}\,dA = -1.564265$  (no fit)")
    ax.axhline(AZ["k2"], color="#c0392b", lw=1.0, ls=":",
               label=r"fitted $k_2=-1.5628$")

    ax.plot(N, keff(N, E), "o", ms=5.5, color="#22313f", zorder=5,
            label="this work: my own minima")
    ax.plot([LN[0]], keff(*LN), "*", ms=17, color="#e08a1e",
            markeredgecolor="#7a4a00", zorder=6,
            label="Lavrov–Nikonov 2026, $N=10^5$ (31 CPU-hours)")

    ax.set_xscale("log")
    ax.set_xlabel("N")
    ax.set_ylabel(r"$k_2^{\rm eff}(N)=-\left[\frac{\pi}{4}N^2-E(N)\right]/N^{3/2}$")
    ax.set_title("The second coefficient is not a fitting parameter", fontsize=12)
    ax.grid(alpha=.25, which="both", lw=.5)
    ax.legend(loc="lower left", fontsize=8.6, framealpha=.95)
    ax.set_ylim(-1.60, -1.30)
    fig.tight_layout()
    fig.savefig(out)
    print("wrote", out)

    print("\n  running coefficient:")
    for n, e in zip(N, E):
        print(f"     N={int(n):5d}   k2_eff = {keff(n,e):.6f}")
    print(f"     N={LN[0]:5d}   k2_eff = {keff(*LN):.6f}   (Lavrov & Nikonov)")
    print(f"     theory                 {K2_THEORY:.6f}")


if __name__ == "__main__":
    import sys
    main(*sys.argv[1:])
