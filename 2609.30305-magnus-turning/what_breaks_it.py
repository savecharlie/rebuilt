"""Eq. (15) is precise. What it is not is robust.

can_you_measure_mu.py found the two-point estimator good to 0.1-2 % at 3 mm
positional noise, which was the opposite of what I expected. So the limit on
this measurement is not variance. This file measures the BIAS from the three
model assumptions the estimator quietly makes, on NOISELESS trajectories, so
that every departure from the true mu is systematic:

    1. the spin axis is exactly vertical
    2. there is no wind
    3. C_L does not depend on the spin parameter S = r omega / |v|

Iris, fire 308.
"""
import numpy as np
from scipy.integrate import solve_ivp
from flight import coeffs, G, BALLS
from can_you_measure_mu import E1, SCENARIOS


def rhs_gen(t, y, kD, kL, what, wind, CLfun, r, omega):
    """Generalised right-hand side: arbitrary spin axis, steady wind, and an
    optional C_L(S). y = [r(3), v(3), psi, s]."""
    v = y[3:6]
    va = v - wind                     # velocity relative to the air
    Va = np.linalg.norm(va)
    if CLfun is not None and Va > 0:
        S = r*omega/Va
        kL = kL*CLfun(S)              # kL passed in as kL/CL_ref
    mg = np.cross(what, va)
    a = np.array([0.0, 0.0, -G]) - kD*Va*va + kL*Va*mg
    W2 = v[0]**2 + v[1]**2
    psidot = (v[0]*a[1] - v[1]*a[0])/W2 if W2 > 0 else 0.0
    return [*v, *a, psidot, np.linalg.norm(v)]


def track_gen(v0, theta, kD, kL, T, f, z0, what, wind, CLfun=None, r=0.0, omega=0.0):
    n = int(round(T*f)) + 1
    ts = np.arange(n)/f
    y0 = [0, 0, z0, v0*np.cos(theta), 0.0, v0*np.sin(theta), 0.0, 0.0]
    s = solve_ivp(rhs_gen, (0, ts[-1]*1.0001), y0,
                  args=(kD, kL, what, wind, CLfun, r, omega),
                  method="DOP853", rtol=1e-12, atol=1e-13, dense_output=True)
    Y = s.sol(ts)
    return ts, Y[:3].T, Y


def axis(eps_deg, phi_deg=0.0):
    e, p = np.radians(eps_deg), np.radians(phi_deg)
    return np.array([np.sin(e)*np.cos(p), np.sin(e)*np.sin(p), np.cos(e)])


if __name__ == "__main__":
    for name in ["fly ball", "free kick"]:
        ball, v0, thd, T, f, W = SCENARIOS[name]
        m, rr, CD, CL = BALLS[ball]
        kD, kL = coeffs(m, rr, CD, CL); mu = CD/CL
        th = np.radians(thd)
        print(f"\n=== {name} ({ball}), true mu = {mu:.3f} "
              f"[noiseless; every deviation below is bias] ===")

        ts, P, Y = track_gen(v0, th, kD, kL, T, f, 2.0, np.array([0, 0, 1.0]),
                             np.zeros(3))
        base = E1(ts, P, W)
        print(f"  control (vertical axis, still air)        mu-hat = {base:.6f}"
              f"   err {abs(base/mu-1):.1e}")

        print("  spin axis tilted from vertical:")
        for eps in [1, 2, 5, 10]:
            worst = 0.0
            for phi in [0, 90, 180, 270]:
                ts, P, _ = track_gen(v0, th, kD, kL, T, f, 2.0,
                                     axis(eps, phi), np.zeros(3))
                e = E1(ts, P, W)
                worst = max(worst, abs(e/mu - 1))
            print(f"    {eps:>3} deg   worst over azimuth: {worst:7.2%} error in mu")

        print("  steady horizontal wind:")
        for u in [0.5, 1.0, 2.0, 5.0]:
            worst = 0.0; arg = None
            for bear in [0, 90, 180, 270]:
                wnd = u*np.array([np.cos(np.radians(bear)), np.sin(np.radians(bear)), 0])
                ts, P, _ = track_gen(v0, th, kD, kL, T, f, 2.0,
                                     np.array([0, 0, 1.0]), wnd)
                e = E1(ts, P, W)
                if abs(e/mu-1) > worst: worst, arg = abs(e/mu-1), (bear, e)
            print(f"    {u:>4.1f} m/s  worst at bearing {arg[0]:>3} deg: "
                  f"mu-hat = {arg[1]:7.4f}  ({worst:7.2%} error)")

        print("  C_L saturating in the spin parameter, C_L = C_L0 /(2+1/S) * 3:")
        for om in [50, 150, 400]:
            # normalise so that the reference spin reproduces the control C_L
            Sref = rr*om/v0
            norm = 1.0/(2 + 1/Sref)
            ts, P, _ = track_gen(v0, th, kD, kL/1.0, T, f, 2.0,
                                 np.array([0, 0, 1.0]), np.zeros(3),
                                 CLfun=lambda S: (1/(2+1/S))/norm, r=rr, omega=om)
            e = E1(ts, P, W)
            print(f"    omega = {om:>4} rad/s (S0 = {Sref:.3f})   "
                  f"mu-hat = {e:7.4f}  ({abs(e/mu-1):7.2%} error)")
