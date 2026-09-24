"""Independent check of the paper's STRUCTURAL claim, not just its energy.

Lavrov & Nikonov report, for N = 100000: sixfold-coordination fraction > 0.86,
<|psi6|> ~ 0.914, and a bulk broken into grains by chains of 5-7 defects.
None of that follows from the energy expansion.  This rebuilds the same
diagnostics on my own configuration, at whatever N I can actually reach, and
draws the orientation map.

Controls, because a bond-order parameter is easy to get silently wrong:
  * a perfect triangular lattice must give |psi6| = 1 exactly in its interior
  * a square lattice must give |psi6| ~ 0
  * Euler: on a Voronoi tessellation of a disk, sum over interior cells of
    (6 - z) must equal 6 minus the boundary correction -- checked as the
    total topological charge of the interior.
"""
import sys
import numpy as np
from scipy.spatial import Voronoi
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def neighbours(xy):
    vor = Voronoi(xy)
    n = len(xy)
    nb = [[] for _ in range(n)]
    for a, b in vor.ridge_points:
        nb[a].append(b)
        nb[b].append(a)
    return nb, vor


def psi6(xy, nb):
    out = np.zeros(len(xy), complex)
    z = np.zeros(len(xy), int)
    for i, js in enumerate(nb):
        if not js:
            continue
        d = xy[js] - xy[i]
        th = np.arctan2(d[:, 1], d[:, 0])
        out[i] = np.exp(6j * th).mean()
        z[i] = len(js)
    return out, z


def control():
    print("=" * 70)
    print("CONTROL — psi6 on lattices whose answer is known")
    print("=" * 70)
    for kind, want in (("triangular", 1.0), ("square", 0.0)):
        k = 9
        i, j = np.meshgrid(np.arange(-k, k + 1), np.arange(-k, k + 1), indexing="ij")
        if kind == "triangular":
            x = (i + 0.5 * j).ravel().astype(float)
            y = (np.sqrt(3) / 2 * j).ravel().astype(float)
        else:
            x, y = i.ravel().astype(float), j.ravel().astype(float)
        xy = np.stack([x, y], 1)
        nb, _ = neighbours(xy)
        p, z = psi6(xy, nb)
        r = np.hypot(x, y)
        core = r < k - 2
        print(f"  {kind:11s}  <|psi6|> interior = {np.abs(p[core]).mean():.6f}"
              f"   (want {want:.1f})   mean z = {z[core].mean():.3f}")


def radial_cdf(xy):
    """Does the configuration actually adopt the arcsine measure?

    My first version binned the radial density and compared bin by bin.  It
    reported median 17 % error and looked like a failure -- but the residuals
    OSCILLATED, which is the signature of discrete concentric rings, not of a
    wrong measure.  The binned density was the wrong instrument.  The
    cumulative distribution averages over the rings:

        F_arcsine(r) = 1 - sqrt(1 - r^2)

    and the scale a sample of N points can resolve at all is 1/sqrt(N).
    """
    r = np.sort(np.hypot(xy[:, 0], xy[:, 1]))
    N = len(r)
    emp = (np.arange(1, N + 1) - 0.5) / N
    arc = 1.0 - np.sqrt(np.clip(1 - r ** 2, 0, 1))
    on_rim = int((r > 1 - 1e-9).sum())
    interior = r < 1 - 1e-9
    return r, emp, arc, on_rim, interior


def analyse(xy, tag):
    nb, vor = neighbours(xy)
    p, z = psi6(xy, nb)
    r = np.hypot(xy[:, 0], xy[:, 1])
    # "bulk" = away from the rim layer of width ~N^{-2/3}... use a generous 0.85
    bulk = r < 0.85
    six = (z[bulk] == 6).mean()
    print(f"\n  {tag}:  N = {len(xy)}")
    print(f"     sixfold fraction in the bulk (r<0.85) = {six:.4f}")
    print(f"     <|psi6|> over all charges             = {np.abs(p).mean():.4f}")
    print(f"     <|psi6|> in the bulk                  = {np.abs(p[bulk]).mean():.4f}")
    print(f"     |<psi6>| (global, phase-coherent)     = {abs(p.mean()):.4f}")
    print(f"     total topological charge sum(6-z) in bulk = {int((6-z[bulk]).sum())}")
    print(f"     Lavrov & Nikonov at N=1e5: >0.86, ~0.914, |<psi6>|~0.095")

    r, emp, arc, on_rim, interior = radial_cdf(xy)
    n = len(r)
    print(f"\n     charges sitting exactly on the rim: {on_rim} of {n} "
          f"({100*on_rim/n:.1f} %)   [2.843 N^(-1/3) = {100*2.84328*n**(-1/3):.1f} %]")
    print(f"     sup|F_emp - F_arcsine| over the interior = "
          f"{np.abs(emp[interior]-arc[interior]).max():.5f}")
    print(f"     1/sqrt(N) = {1/np.sqrt(n):.5f}  <- what a sample this size can resolve")
    print("        r      empirical F    arcsine F      diff")
    for q in (0.1, 0.25, 0.5, 0.7):
        i = int(q * n)
        print(f"     {r[i]:6.4f}   {emp[i]:11.5f}  {arc[i]:11.5f}   {emp[i]-arc[i]:+9.5f}")
    print("     (the deficit is the mass the continuum puts near the rim and the")
    print("      discrete system puts ON it; no density can represent that layer.)")

    fig, axes = plt.subplots(1, 2, figsize=(11.6, 5.8), dpi=170)
    ang = (np.angle(p) % (2 * np.pi)) / 6.0
    s = max(1.0, 2600.0 / len(xy))
    a0 = axes[0].scatter(xy[:, 0], xy[:, 1], c=np.degrees(ang), s=s,
                         cmap="hsv", vmin=0, vmax=60)
    axes[0].set_title(r"local lattice orientation  $\arg\psi_6/6$  (deg mod 60)",
                      fontsize=10)
    plt.colorbar(a0, ax=axes[0], fraction=.046)
    cols = np.where(z == 5, "#c0392b", np.where(z == 7, "#2471a3", "#d8d8d8"))
    axes[1].scatter(xy[:, 0], xy[:, 1], c=cols, s=s * 1.3)
    axes[1].set_title("five-fold (red) and seven-fold (blue) charges", fontsize=10)
    for ax in axes:
        ax.set_aspect("equal")
        ax.add_patch(plt.Circle((0, 0), 1, fill=False, lw=.8, color="k"))
        ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle(f"a disk cannot be one crystal — N = {len(xy)}, my own minimum",
                 fontsize=11)
    fig.tight_layout()
    out = f"structure_N{len(xy)}.png"
    fig.savefig(out)
    print(f"     wrote {out}")


def main(N=1350):
    control()
    from ladder import best
    N = int(N)
    print(f"\n  optimising N = {N} for the structure map ...", flush=True)
    E, xy, Nb = best(N, seeds=1, maxiter=1200)
    np.save(f"config_N{N}.npy", xy)
    print(f"     E = {E:.8f}   Nb = {Nb}")
    analyse(xy, "my minimum")


if __name__ == "__main__":
    main(*sys.argv[1:])
