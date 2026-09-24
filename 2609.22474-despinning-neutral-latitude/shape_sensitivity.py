"""What the neutral latitude actually depends on: the degree-4 part of the figure.

structure_robust.py shows the area-neutral latitude is pinned near the magic
angle for a spheroid and collapses for a point-mass (Roche) figure.  The two
differ by 21 degrees at a/c = 1.49.  This file finds the knob.

Write any axisymmetric figure as

    r(phi) / R = 1 + s2 P2(sin phi) + s4 P4(sin phi) + ...

A rotational perturbation is pure degree 2 at first order, so s2 carries the
flattening and s4 measures how non-spheroidal the figure is.  Measured here:

  1. CONTROL: the numerical pipeline (np.gradient on a sampled shape) must
     reproduce the analytic closed forms of neutral.py.
  2. the area-neutral latitude as a function of s4/s2 at fixed s2
  3. the s4/s2 that the two exact figures actually carry

Run:  python3 shape_sensitivity.py
"""
import sys, pathlib
import numpy as np
from scipy import integrate, optimize, special

sys.path.insert(0, str(pathlib.Path(__file__).parent))
sys.path.insert(0, str(pathlib.Path(__file__).parents[1] / "tools"))
from figures import spheroid_radius, spheroid_axes, line_elements
from neutral import phi_lon_instant, phi_lat_instant
from structure_robust import PHI, roche_at_ratio, spheroid_at_ratio, neutral, normalise, DEG
import dull

X = np.sin(PHI)
P = {n: special.eval_legendre(n, X) for n in (0, 2, 4, 6, 8)}


def legendre_coeffs(r, nmax=8):
    """s_n of r(phi)/R = sum s_n P_n(sin phi), by orthogonality over mu=sin(phi).
    The integrand is even in mu, so integrating 0..1 and using the full-interval
    norm 2/(2n+1) is correct for even n only -- which is all a figure has."""
    mu = np.sin(PHI)
    o = np.argsort(mu)
    out = {}
    for n in (0, 2, 4, 6, 8):
        Pn = special.eval_legendre(n, mu)
        out[n] = (2 * n + 1) * integrate.simpson((r * Pn)[o], x=mu[o])
    return {n: out[n] / out[0] for n in (2, 4, 6, 8)}


def figure(s2, s4):
    r = 1.0 + s2 * P[2] + s4 * P[4]
    return normalise(r)


if __name__ == "__main__":
    print(__doc__.split("Run:")[0])
    with dull.watch() as w:
        # ------------------------------------------------- 1. control
        print("=" * 76)
        print("1. CONTROL: numerical pipeline vs the analytic closed forms")
        # The brute-force "instantaneous" strain is a one-sided finite difference in
        # a/c, so it carries an O(step) offset.  dull.scaling() measures that
        # exponent as 0.9995 over four decades of step, i.e. it is the step and not
        # a fault; step 1e-5 puts the offset at ~1e-4 deg.
        DD = 1.00001
        worst = 0.0
        for ratio in (1.05, 1.2, 1.4):
            f = 1 - 1 / ratio
            n = neutral(spheroid_at_ratio(ratio), spheroid_at_ratio(ratio * DD))
            e = (float(phi_lon_instant(f)), float(phi_lat_instant(f)))
            worst = max(worst, abs(n[0] - e[0]), abs(n[1] - e[1]))
            print(f"   a/c={ratio:4.2f}  lon num {n[0]:8.4f} exact {e[0]:8.4f} |"
                  f"  lat num {n[1]:8.4f} exact {e[1]:8.4f}")
        print(f"   worst {worst:.2e} deg")
        assert worst < 1e-3, "numerical pipeline disagrees with the closed form"

        # ------------------------------------------------- 2. sensitivity to s4
        print()
        print("=" * 76)
        print("2. area-neutral latitude vs s4/s2, at fixed s2 (a/c held near 1.3)")
        a, c = spheroid_axes(1 - 1 / 1.3, R=1.0)
        s_sph = legendre_coeffs(spheroid_radius(PHI, a, c))
        s2 = s_sph[2]
        print(f"   spheroid at a/c=1.30 has s2 = {s2:+.5f}, s4/s2 = {s_sph[4]/s2:+.5f}")
        print(f"{'s4/s2':>9} {'lon':>9} {'lat':>9} {'area':>9}")
        ratios, areas = [], []
        for k in (-0.30, -0.20, -0.10, 0.0, 0.10, 0.20, 0.30):
            r1 = figure(s2, k * s2)
            r2 = figure(s2 * 1.0005, k * s2 * 1.0005)
            n = neutral(r1, r2)
            ratios.append(k); areas.append(n[2])
            print(f"{k:+9.2f} {n[0]:9.4f} {n[1]:9.4f} {n[2]:9.4f}")
        slope = np.polyfit(ratios, areas, 1)[0]
        print(f"   d(area-neutral)/d(s4/s2) = {slope:+.2f} deg per unit")

        # ------------------------------------------------- 3. what the real figures carry
        print()
        print("=" * 76)
        print("3. s4/s2 of the two exact figures, and the shift it predicts")
        print(f"{'a/c':>6} | {'spheroid s4/s2':>15} {'Roche s4/s2':>13} |"
              f" {'predicted d(phi)':>17} {'measured d(phi)':>16}")
        for ratio in (1.1, 1.2, 1.3, 1.4, 1.49):
            asp, csp = spheroid_axes(1 - 1 / ratio, R=1.0)
            ssp = legendre_coeffs(spheroid_radius(PHI, asp, csp))
            _, rr = roche_at_ratio(ratio)
            sro = legendre_coeffs(rr)
            ksp, kro = ssp[4] / ssp[2], sro[4] / sro[2]
            nsp = neutral(spheroid_at_ratio(ratio), spheroid_at_ratio(ratio * 1.0005))
            _, r2 = roche_at_ratio(ratio * 1.0005)
            nro = neutral(rr, r2)
            print(f"{ratio:6.2f} | {ksp:+15.5f} {kro:+13.5f} |"
                  f" {slope*(kro-ksp):+17.3f} {nro[2]-nsp[2]:+16.3f}")
    print(f"\ndull.watch: {len(w.findings)} FP events")
    for ev in w.findings[:6]:
        print("   ", ev)
