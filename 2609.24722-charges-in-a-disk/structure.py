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


def radial_density(xy, nbins=26):
    """Does the configuration actually adopt the arcsine measure?

    Everything in k1 and k2 rests on rho_eq(r) = 1/(2 pi sqrt(1-r^2)).  That is
    a claim about where the charges go, and it can be read straight off a
    configuration instead of assumed.
    """
    r = np.hypot(xy[:, 0], xy[:, 1])
    N = len(xy)
    edges = np.linspace(0, 1, nbins + 1)
    cnt, _ = np.histogram(r, bins=edges)
    area = np.pi * (edges[1:] ** 2 - edges[:-1] ** 2)
    meas = cnt / area / N
    mid = 0.5 * (edges[1:] + edges[:-1])
    # predicted: average of rho_eq over the annulus = [sqrt(1-a^2)-sqrt(1-b^2)]/area
    pred = (np.sqrt(np.clip(1 - edges[:-1] ** 2, 0, 1))
            - np.sqrt(np.clip(1 - edges[1:] ** 2, 0, 1))) / area
    return mid, meas, pred


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

    mid, meas, pred = radial_density(xy)
    ok = mid < 0.93
    rel = np.abs(meas[ok] - pred[ok]) / pred[ok]
    print(f"     radial density vs arcsine measure (r<0.93): "
          f"median |rel| = {np.median(rel):.4f}, max = {rel.max():.4f}")
    print("        r    measured    arcsine     ratio")
    for m, a, b in list(zip(mid, meas, pred))[::4]:
        print(f"     {m:5.3f}  {a:9.4f}  {b:9.4f}   {a/b:7.4f}")

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
