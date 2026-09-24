"""The neutral latitude of despinning strain: exact, and checked against brute force.

Lock (2026) computes, for each of a sequence of HERCULES equilibrium figures, the
change in the local line elements (his Eqs. 3-5) at fixed latitude, and reports
that the sign of the area change flips "poleward of about 20-40 degrees".

That band has a closed form.  For a volume-preserving oblate spheroid of
flattening f (u = 1 - f = c/a):

  LONGITUDINAL (E-W, his Eq. 4):     tan^2 phi_0 = u^2 / 2
  LATITUDINAL  (N-S, his Eq. 3):     tan^2 phi_0 = u^2 [7 - 8u^2 + sqrt(64u^4 - 104u^2 + 49)] / 4

Both reduce to tan^2 phi = 1/2, phi = 35.2644 deg, at u = 1 -- the zero of the
Legendre polynomial P2, which is where they have to be, because a rotational
figure perturbation is pure degree 2.  Flattening splits them, and it splits them
in OPPOSITE directions.

For a FINITE change from u1 to u2 (which is what his Fig. 4 plots -- cumulative,
not instantaneous) the longitudinal condition is still exactly solvable:

  tan^2 phi_0 = 1 / [A B (A + B)],   A = u1^{-2/3}, B = u2^{-2/3}

which collapses to u^2/2 when u1 = u2.

Run:  python3 neutral.py
"""
import sys, pathlib
import numpy as np
from scipy import optimize

sys.path.insert(0, str(pathlib.Path(__file__).parent))
sys.path.insert(0, str(pathlib.Path(__file__).parents[1] / "tools"))
from figures import spheroid_radius, spheroid_drdphi, spheroid_axes, line_elements
import dull

DEG = 180.0 / np.pi

# ------------------------------------------------------------------ closed forms

def phi_lon_instant(f):
    """Instantaneous longitudinal neutral latitude, degrees."""
    u = 1.0 - np.asarray(f, float)
    return np.arctan(u / np.sqrt(2.0)) * DEG


def phi_lat_instant(f):
    """Instantaneous latitudinal neutral latitude, degrees."""
    u = 1.0 - np.asarray(f, float)
    u2 = u * u
    disc = 64 * u2 * u2 - 104 * u2 + 49        # discriminant in u^2 is -1728 < 0: always > 0
    return np.arctan(np.sqrt(u2 * (7 - 8 * u2 + np.sqrt(disc)) / 4.0)) * DEG


def phi_lon_cumulative(f1, f2):
    """Neutral latitude of the FINITE longitudinal change from f1 to f2, degrees."""
    A = (1.0 - f1) ** (-2.0 / 3.0)
    B = (1.0 - f2) ** (-2.0 / 3.0)
    return np.arctan(np.sqrt(1.0 / (A * B * (A + B)))) * DEG

# ------------------------------------------------------- brute force, the paper's way

def strains(f1, f2, n=4001):
    """Fractional change in dl_lat, dl_lon and dA from flattening f1 to f2,
    computed exactly as Lock's Eqs. 3-5, at fixed geocentric latitude."""
    phi = np.linspace(1e-9, np.pi / 2 - 1e-9, n)
    out = []
    for f in (f1, f2):
        a, c = spheroid_axes(f, R=1.0)
        r = spheroid_radius(phi, a, c)
        out.append(line_elements(phi, r, spheroid_drdphi(phi, a, c)))
    (lat1, lon1), (lat2, lon2) = out
    return phi, lat2 / lat1 - 1.0, lon2 / lon1 - 1.0, (lat2 * lon2) / (lat1 * lon1) - 1.0


def _zero_crossing(phi, y):
    i = np.where(np.sign(y[:-1]) != np.sign(y[1:]))[0]
    if len(i) != 1:
        return np.nan
    i = i[0]
    t = -y[i] / (y[i + 1] - y[i])
    return (phi[i] + t * (phi[i + 1] - phi[i])) * DEG


if __name__ == "__main__":
    print(__doc__.split("Run:")[0])
    with dull.watch() as w:
        print("=" * 78)
        print("1. INSTANTANEOUS closed forms vs brute force (df = 1e-6)")
        print(f"{'f':>6} {'a/c':>7} | {'lon exact':>10} {'lon brute':>10} | {'lat exact':>10} {'lat brute':>10}")
        worst = 0.0
        for f in (0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6):
            phi, elat, elon, _ = strains(f, f + 1e-6)
            bl, bt = _zero_crossing(phi, elon), _zero_crossing(phi, elat)
            xl, xt = float(phi_lon_instant(f)), float(phi_lat_instant(f))
            worst = max(worst, abs(bl - xl), abs(bt - xt))
            print(f"{f:6.2f} {1/(1-f):7.4f} | {xl:10.5f} {bl:10.5f} | {xt:10.5f} {bt:10.5f}")
        print(f"   worst disagreement: {worst:.2e} deg")
        assert worst < 2e-3, "closed form disagrees with brute force"

        print()
        print("=" * 78)
        print("2. The two neutral latitudes COINCIDE at f=0 and SPLIT in opposite directions")
        print(f"   f=0     : lon {phi_lon_instant(0.0):.4f}  lat {phi_lat_instant(0.0):.4f}"
              f"   (arctan(1/sqrt2) = {np.arctan(1/np.sqrt(2))*DEG:.6f}, the zero of P2)")
        print(f"   f=0.5   : lon {phi_lon_instant(0.5):.4f}  lat {phi_lat_instant(0.5):.4f}")
        print(f"   arcsin(1/3) = {np.arcsin(1/3)*DEG:.6f} deg  <- the f=0.5 longitudinal value, exactly")
        assert abs(float(phi_lon_instant(0.5)) - np.arcsin(1 / 3) * DEG) < 1e-9

        m = optimize.minimize_scalar(lambda f: -float(phi_lat_instant(f)), bounds=(0.0, 0.9),
                                     method="bounded", options={"xatol": 1e-10})
        print(f"   latitudinal neutral latitude PEAKS at f = {m.x:.6f} (a/c = {1/(1-m.x):.4f}),"
              f" phi = {-m.fun:.4f} deg")

        print()
        print("=" * 78)
        print("3. CUMULATIVE change to the present day (f_now = 1/298.257)")
        f_now = 1 / 298.257
        print(f"{'f_init':>7} {'a/c':>7} | {'lon exact':>10} {'lon brute':>10} | {'lat brute':>10} {'area brute':>10}")
        for f0 in (0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5):
            phi, elat, elon, earea = strains(f0, f_now)
            print(f"{f0:7.3f} {1/(1-f0):7.4f} | {float(phi_lon_cumulative(f0,f_now)):10.5f}"
                  f" {_zero_crossing(phi,elon):10.5f} | {_zero_crossing(phi,elat):10.5f}"
                  f" {_zero_crossing(phi,earea):10.5f}")

        print()
        print("=" * 78)
        print("4. Magnitudes: cumulative fractional change in local AREA")
        print(f"{'f_init':>7} | {'equator':>9} {'45 deg':>9} {'pole':>9}")
        for f0 in (0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5):
            phi, _, _, earea = strains(f0, f_now)
            g = lambda d: earea[np.argmin(abs(phi * DEG - d))]
            print(f"{f0:7.3f} | {g(0)*100:8.2f}% {g(45)*100:8.2f}% {g(90)*100:8.2f}%")
    print(f"\ndull.watch: {len(w.findings)} FP events")
    for ev in w.findings[:5]:
        print("   ", ev)
