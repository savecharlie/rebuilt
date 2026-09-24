"""Validate the instrument BEFORE it says anything new.

Every check here has an answer known independently of this rebuild: a closed
form, a textbook constant, or a measured property of the Earth.  If any of
these is red, nothing downstream is worth reading.

Run:  python3 validate.py
"""
import sys, pathlib
import numpy as np
from scipy import integrate, optimize

sys.path.insert(0, str(pathlib.Path(__file__).parent))
sys.path.insert(0, str(pathlib.Path(__file__).parents[1] / "tools"))
from figures import (spheroid_radius, spheroid_drdphi, spheroid_axes, line_elements,
                     spheroid_volume, spheroid_area, meridian_quarter_arc,
                     maclaurin_omega2, roche_surface, radau_darwin_f)
import dull

RESULTS = []

def check(name, got, want, tol, note=""):
    err = abs(got - want) / (abs(want) if want else 1.0)
    ok = err <= tol
    RESULTS.append((ok, name, got, want, err, note))
    print(f"[{'ok ' if ok else 'FAIL'}] {name:<52s} got {got:< 18.10g} want {want:< 18.10g} rel {err:.2e}  {note}")
    return ok

# ---------------------------------------------------------------- 1. geometry engine
# The line elements of Eq. 3/4 are integrated and must reproduce the closed-form
# area and meridian arc of the spheroid.  This is the whole engine; if it is
# wrong every strain below is wrong.
a, c = 1.3, 0.8
phi = np.linspace(0, np.pi / 2, 200001)
r = spheroid_radius(phi, a, c)
dr = spheroid_drdphi(phi, a, c)
dl_lat, dl_lon = line_elements(phi, r, dr)

check("meridian quarter-arc from Eq.3", integrate.simpson(dl_lat, x=phi),
      meridian_quarter_arc(a, c), 1e-10, "vs a*E(e^2)")
check("surface area from Eqs.3-5", 4 * np.pi * integrate.simpson(dl_lat * dl_lon, x=phi),
      spheroid_area(a, c), 1e-10, "vs exact oblate area")

# analytic dr/dphi against finite difference -- catches a sign or factor slip
h = 1e-6
fd = (spheroid_radius(phi[1:-1] + h, a, c) - spheroid_radius(phi[1:-1] - h, a, c)) / (2 * h)
check("analytic dr/dphi", float(np.max(np.abs(fd - dr[1:-1]))), 0.0, 1.0, "max abs diff (see next)")
assert np.max(np.abs(fd - dr[1:-1])) < 1e-7, "analytic derivative disagrees with FD"

# volume-preserving axes
for f in (0.0, 0.1, 0.3, 0.5):
    aa, cc = spheroid_axes(f, R=1.0)
    check(f"volume preserved at f={f}", spheroid_volume(aa, cc), 4 / 3 * np.pi, 1e-12,
          f"a/c={aa/cc:.4f}")

# ---------------------------------------------------------------- 2. Maclaurin sequence
# Three textbook constants of the homogeneous rotating sequence.
res = optimize.minimize_scalar(lambda e: -maclaurin_omega2(e), bounds=(0.5, 0.999), method="bounded",
                               options={"xatol": 1e-12})
check("Maclaurin: e at max omega^2", -res.x * -1, 0.9299514, 1e-5, "Chandrasekhar")
check("Maclaurin: max omega^2/(pi G rho)", -res.fun, 0.4493, 1e-3, "")
# small-e limit must be f = 5q/4 with q = omega^2 a^3/GM = omega^2/((4/3)pi G rho)
# NOTE e=1e-4 makes this check read inf: maclaurin_omega2 is a difference of two
# terms of size 6e8 whose difference is 5e-9, so float64 cancels it to zero.  A
# dull fault in the validation itself, found by the validation itself.  e=0.01
# leaves ~7 digits, which is all this check needs.
e = 1e-2
f_small = 1 - np.sqrt(1 - e ** 2)
q_small = maclaurin_omega2(e) * 3.0 / 4.0        # omega^2/(pi G rho) -> omega^2 a^3/GM
check("Maclaurin small-e: f/q", f_small / q_small, 1.25, 1e-3, "= 5/4 exactly")

# ---------------------------------------------------------------- 3. Darwin-Radau
q_E = (7.292115e-5) ** 2 * 6378137.0 ** 3 / 3.986004418e14
check("Darwin-Radau homogeneous lam", radau_darwin_f(q_E, 0.4) / q_E, 1.25, 1e-12, "")
f_E = radau_darwin_f(q_E, 0.3307)
check("Earth f from C/MR^2=0.3307", 1 / f_E, 298.257, 5e-3, "observed 1/298.257 (WGS84)")

# ---------------------------------------------------------------- 4. Roche model
# Point mass + rotation.  Two exact limits: q->0 gives a sphere, and the
# equatorial radius satisfies the cubic exactly at the break-up value q=8/27
# where r_eq/r_pol = 3/2 (the classical Roche limit of the model).
check("Roche: q->0 is a sphere", float(roche_surface(np.array([0.0]), 1e-12)[0]), 1.0, 1e-9, "")
check("Roche: break-up q=8/27 -> a/c = 3/2", float(roche_surface(np.array([0.0]), 8 / 27)[0]),
      1.5, 1e-7, "double root: np.roots loses half the digits")

# ---------------------------------------------------------------- 5. dull.py is live
with dull.watch() as w:
    _ = spheroid_radius(np.linspace(0, np.pi / 2, 1001), *spheroid_axes(0.5))
print(f"\ndull.watch: {len(w.findings)} FP events during the geometry pass")

nbad = sum(1 for ok, *_ in RESULTS if not ok)
print(f"\n{len(RESULTS) - nbad}/{len(RESULTS)} validation checks pass")
sys.exit(1 if nbad else 0)
