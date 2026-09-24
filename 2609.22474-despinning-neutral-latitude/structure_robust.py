"""Does the neutral latitude depend on the interior structure?

The closed forms in neutral.py are for an oblate SPHEROID -- exactly the figure
of a homogeneous (Maclaurin) body.  Earth is centrally condensed (C/MR^2 = 0.331
against 0.4 for uniform), so its equipotential surface is not a spheroid: it
carries degree-4 and higher terms, and HERCULES resolves them (Lock truncates at
degree 12).  So the question the closed form cannot answer by itself is whether
the neutral latitude is a property of the geometry or of the interior.

Two exact extremes bracket every real body:

  MACLAURIN  uniform density, C/MR^2 = 2/5, surface is exactly a spheroid
  ROCHE      all mass at the centre, C/MR^2 -> 0, surface from GM/r + w^2 r^2 cos^2/2

Compared here at EQUAL a/c and equal enclosed volume.  The Roche model cannot
reach a/c > 3/2 (it breaks up at q = 8/27), which is itself worth saying out loud
about a paper whose extreme case is a/c = 2.

Run:  python3 structure_robust.py
"""
import sys, pathlib
import numpy as np
from scipy import integrate, optimize

sys.path.insert(0, str(pathlib.Path(__file__).parent))
sys.path.insert(0, str(pathlib.Path(__file__).parents[1] / "tools"))
from figures import spheroid_radius, spheroid_drdphi, spheroid_axes, line_elements
import dull

DEG = 180.0 / np.pi
PHI = np.linspace(1e-7, np.pi / 2 - 1e-7, 8001)


def roche_shape(q):
    """r(phi)/R_pol for the Roche figure, from the exact cubic, on the PHI grid."""
    c2 = np.cos(PHI) ** 2
    x = np.empty_like(PHI)
    for i, cc in enumerate(c2):
        if q * cc < 1e-14:
            x[i] = 1.0
            continue
        # q cc x^3/2 - x + 1 = 0, small positive root, solved by Newton from x=1
        g = lambda z: q * cc * z ** 3 / 2 - z + 1
        dg = lambda z: 1.5 * q * cc * z ** 2 - 1
        x[i] = optimize.brentq(g, 1.0, 1.5000000001, xtol=1e-15, rtol=8.9e-16)
    return x


def volume(r):
    return 4 * np.pi / 3 * integrate.simpson(r ** 3 * np.cos(PHI), x=PHI)


def normalise(r):
    return r * (4 * np.pi / 3 / volume(r)) ** (1 / 3)


def elements(r):
    dr = np.gradient(r, PHI, edge_order=2)
    return line_elements(PHI, r, dr)


def roche_at_ratio(ratio):
    """q such that a/c = ratio; returns the volume-normalised shape."""
    def err(q):
        x = roche_shape(q)
        return x[0] / x[-1] - ratio          # PHI[0] ~ equator, PHI[-1] ~ pole
    q = optimize.brentq(err, 1e-9, 8 / 27 - 1e-12, xtol=1e-14)
    return q, normalise(roche_shape(q))


def spheroid_at_ratio(ratio):
    a, c = spheroid_axes(1 - 1 / ratio, R=1.0)
    return spheroid_radius(PHI, a, c)


def _cross(y):
    i = np.where(np.sign(y[:-1]) != np.sign(y[1:]))[0]
    if len(i) != 1:
        return np.nan
    i = i[0]
    return (PHI[i] - y[i] * (PHI[i + 1] - PHI[i]) / (y[i + 1] - y[i])) * DEG


def neutral(r1, r2):
    (lat1, lon1), (lat2, lon2) = elements(r1), elements(r2)
    return (_cross(lon2 / lon1 - 1), _cross(lat2 / lat1 - 1),
            _cross(lat2 * lon2 / (lat1 * lon1) - 1))


if __name__ == "__main__":
    print(__doc__.split("Run:")[0])
    with dull.watch() as w:
        d = 1.0005
        print(f"{'a/c':>6} | {'MACLAURIN (spheroid)':^34} | {'ROCHE (point mass)':^34}")
        print(f"{'':>6} | {'lon':>10} {'lat':>10} {'area':>10} | {'lon':>10} {'lat':>10} {'area':>10}")
        rows = []
        for ratio in (1.02, 1.05, 1.1, 1.2, 1.3, 1.4, 1.49):
            s1, s2 = spheroid_at_ratio(ratio), spheroid_at_ratio(ratio * d)
            m = neutral(s1, s2)
            _, r1 = roche_at_ratio(ratio)
            _, r2 = roche_at_ratio(ratio * d)
            rr = neutral(r1, r2)
            rows.append((ratio, m, rr))
            print(f"{ratio:6.2f} | {m[0]:10.4f} {m[1]:10.4f} {m[2]:10.4f} |"
                  f" {rr[0]:10.4f} {rr[1]:10.4f} {rr[2]:10.4f}")

        print()
        print("area-neutral latitude, Maclaurin vs Roche, same a/c:")
        for ratio, m, rr in rows:
            print(f"   a/c={ratio:5.2f}  spheroid {m[2]:8.4f}   Roche {rr[2]:8.4f}"
                  f"   diff {rr[2]-m[2]:+7.4f} deg")

        print()
        print("Roche model ceiling: the figure of a fully centrally-condensed body")
        q = 8 / 27
        x = roche_shape(q - 1e-12)
        print(f"   break-up at q = 8/27 = {q:.8f}, a/c = {x[0]/x[-1]:.8f}  (exactly 3/2)")
        print("   -> a body this condensed CANNOT reach the a/c = 2 of the paper's 2.7 L_EM case.")
    print(f"\ndull.watch: {len(w.findings)} FP events")
    for ev in w.findings[:6]:
        print("   ", ev)
