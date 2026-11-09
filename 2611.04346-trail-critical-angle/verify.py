#!/usr/bin/env python3
"""Rebuild equation 2 of arXiv:2611.04346, and carry it one step the paper does not take.

The paper's equation 2 is

    theta_c = arccos(C_trail / C_ground)

the angle, measured from the trail's own direction, at which a walker should leave a
straight trail and strike out for the target. It is stated, used, and verified against
simulation -- it is not derived in the text, and the word that would place it is absent:
I grepped for snell, refract, fermat, head wave, optic. "Snell" appears once, as
Doyle & Snell, *Random walks and electric networks*; "seismic" once, in a reference
about survey lines used as animal trails. Both are coincidences.

It is the seismic head wave. A trail is a fast refractor embedded in a slow medium
(cost = slowness, so C == refractive index). A ray entering the refractor runs along it
and radiates back out at the critical angle, and the head-wave traveltime is
t = x/V2 + 2z*sqrt(V2^2-V1^2)/(V1*V2) -- the same shape as the cost derived below.
With the angle measured from the NORMAL rather than the tangent, sin(theta_normal) =
V_slow/V_fast = C_trail/C_ground, which is equation 2's complement.

Three checks, cheapest and most certain first:

  1. the closed form for the optimal detour cost, against brute force
  2. theta_c itself, measured off the brute-force optimum, against equation 2
  3. the CROSSOVER DISTANCE -- the quantity the paper has no name for and which is the
     one an empiricist can actually measure in the field: how far apart must two
     points lying a perpendicular distance h off a trail be before the trail is worth
     using at all?

Iris (Opus 5, 1M), fire 416, 9 Nov 2026.
"""
import math

C_T = 1.0                     # Table 1: trail traversal cost is fixed at 1.0
RATIOS = [1.25, 1.5, 2.0, 3.0, 5.0, 8.0, 20.0, 65.0]   # Table 1: C_ground in [1.0, 65.0]


def brute_one_sided(Cg, d, h, n=400001):
    """Start ON the trail at the origin; target at (d, h). Minimise over the leaving
    point x by dense scan. Returns (cost, x*)."""
    best, bx = float("inf"), None
    for i in range(n):
        x = d * i / (n - 1)
        c = C_T * x + Cg * math.hypot(d - x, h)
        if c < best:
            best, bx = c, x
    return best, bx


def closed_one_sided(Cg, d, h):
    """C_T*d + h*sqrt(Cg^2 - C_T^2), valid while the optimum is interior."""
    return C_T * d + h * math.sqrt(Cg * Cg - C_T * C_T)


def brute_two_sided(Cg, d, h, n=1201):
    """Both endpoints a perpendicular distance h off the trail, separated by d along it.
    Minimise over the entry and exit points independently -- NOT assuming the walker
    uses the trail at all (entry beyond exit degenerates to the straight line)."""
    best = Cg * d                                   # the straight line, h0 == h1
    for i in range(n):
        a = d * i / (n - 1)
        for j in range(i, n):
            b = d * j / (n - 1)
            c = Cg * math.hypot(a, h) + C_T * (b - a) + Cg * math.hypot(d - b, h)
            if c < best:
                best = c
    return best


def d_cross(Cg, h):
    """2h*sqrt((r+1)/(r-1)) -- the along-trail separation at which using the trail first
    pays. The seismic crossover distance, in trail coordinates."""
    r = Cg / C_T
    return 2.0 * h * math.sqrt((r + 1.0) / (r - 1.0))


print("1 & 2.  one-sided detour: closed form and the angle it implies")
print("        start on the trail at the origin, target at (d, h) = (40, 10)\n")
print("   Cg/Ct     brute cost   closed form    rel err      theta* measured   eq2     diff")
d, h = 40.0, 10.0
bad = 0
for Cg in RATIOS:
    bc, bx = brute_one_sided(Cg, d, h)
    cc = closed_one_sided(Cg, d, h)
    u = d - bx                                   # along-trail distance still to run
    th = math.degrees(math.atan2(h, u)) if u > 0 else 90.0
    eq2 = math.degrees(math.acos(C_T / Cg))
    rel = abs(bc - cc) / cc
    bad += (rel > 2e-6) or (abs(th - eq2) > 0.05)
    print("  %6.2f   %10.4f   %10.4f   %9.2e      %7.3f deg   %7.3f  %+6.3f"
          % (Cg, bc, cc, rel, th, eq2, th - eq2))

print("\n   (the closed form is C_T*d + h*sqrt(Cg^2 - C_T^2); the measured angle is read")
print("    off the brute-force optimum, never assumed. The interior optimum needs")
print("    d >= h*C_T/sqrt(Cg^2-C_T^2), which is %.3f here at Cg=1.25.)"
      % (h * C_T / math.sqrt(1.25 ** 2 - 1)))

print("\n3.  crossover: both points a distance h = 10 off the trail")
print("    d_cross = 2h*sqrt((r+1)/(r-1)).  Below it the straight line wins and the")
print("    trail captures nothing; above it the trail route is strictly cheaper.\n")
print("   Cg/Ct    d_cross   cost at 0.9*d_cross       cost at 1.1*d_cross")
print("                      trail      straight       trail      straight")
for Cg in RATIOS:
    dc = d_cross(Cg, 10.0)
    row = [dc]
    for f in (0.9, 1.1):
        dd = f * dc
        row += [brute_two_sided(Cg, dd, 10.0), Cg * dd]
    lo_t, lo_s, hi_t, hi_s = row[1:]
    ok = (lo_t >= lo_s - 1e-9) and (hi_t < hi_s - 1e-9)
    bad += not ok
    print("  %6.2f  %8.2f   %9.3f %9.3f     %9.3f %9.3f   %s"
          % (Cg, dc, lo_t, lo_s, hi_t, hi_s, "ok" if ok else "FAIL"))

print("\n  r -> 1   : d_cross -> infinity (a trail no better than the ground captures nothing)")
print("  r -> inf : d_cross -> 2h        (even a free trail needs a separation of twice the offset)")
for r in (1.01, 1.1, 2.0, 10.0, 65.0, 1e6):
    print("    r = %-9g d_cross/h = %.4f" % (r, d_cross(r, 1.0)))

print("\n%s" % ("ALL CHECKS PASS" if bad == 0 else "%d CHECK(S) FAILED" % bad))
raise SystemExit(1 if bad else 0)
