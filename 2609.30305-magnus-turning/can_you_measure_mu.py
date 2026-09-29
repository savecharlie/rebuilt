"""Can Eq. (15) actually be used?

The paper offers

    C_D/C_L = ln(|w_0|/|w_f|) / dpsi                                    (15)

and calls it a measurement: "the lift-to-drag ratio of a real ball can be read
off a single tracked trajectory". The identity is exact -- I confirmed it to
1e-12 in check_invariants.py. This file asks the other question, which the paper
does not: given a real tracker, what is the error bar?

Two estimators on the same synthetic data:

  E1  the paper's two-point estimator. Local quadratic fits to the first and
      last W frames give the horizontal velocity at each end; mu-hat is the log
      ratio over the heading swing. Uses no knowledge of g, m, A, rho, or even
      which ball it is.

  E2  full nonlinear least squares for (r0, v0, kD, kL), eight parameters,
      assuming the model of Eq. (1) exactly and knowing g. mu-hat = kD/kL.
      Its covariance is bounded below by the Cramer-Rao bound, which is also
      computed analytically here.

Iris, fire 308.
"""
import numpy as np
from scipy.optimize import least_squares
from flight import fly, rhs, coeffs, G, BALLS
from scipy.integrate import solve_ivp

RNG = np.random.default_rng(20260929)


def track(v0, theta, kD, kL, T, f, z0=0.0):
    """True sampled positions: N frames at rate f over duration T."""
    n = int(round(T*f)) + 1
    ts = np.arange(n)/f
    sol = fly(v0, theta, kD, kL, z0=z0, t_max=T*1.001, rtol=1e-12, atol=1e-13)
    Y = sol.sol(ts)
    return ts, Y[:3].T, Y          # times, positions (n,3), full state


def E1(ts, P, W):
    """Two-point estimator of mu from noisy positions."""
    out = []
    for sl in (slice(0, W), slice(len(ts)-W, len(ts))):
        t = ts[sl] - ts[sl].mean()
        A = np.c_[np.ones_like(t), t, t*t]
        c, *_ = np.linalg.lstsq(A, P[sl, :2], rcond=None)
        vx, vy = c[1]                       # derivative at the window centre
        out.append((np.hypot(vx, vy), np.arctan2(vy, vx)))
    (wa, pa), (wb, pb) = out
    dpsi = np.unwrap([pa, pb])[1] - pa
    if abs(dpsi) < 1e-12:
        return np.nan
    return np.log(wa/wb)/dpsi


def _model(p, ts):
    r0 = p[:3]; v0 = p[3:6]; kD, kL = np.exp(p[6]), np.exp(p[7])
    y0 = list(r0) + list(v0) + [np.arctan2(v0[1], v0[0]), 0.0]
    s = solve_ivp(rhs, (0, ts[-1]*1.0001), y0, args=(kD, kL, +1),
                  method="DOP853", rtol=1e-11, atol=1e-12, dense_output=True)
    return s.sol(ts)[:3].T


def E2(ts, P, p_guess):
    r = least_squares(lambda p: (_model(p, ts) - P).ravel(), p_guess,
                      xtol=1e-14, ftol=1e-14, gtol=1e-14)
    return np.exp(r.x[6] - r.x[7]), r.x     # mu = kD/kL


def crb_sigma_mu(p_true, ts, sigma):
    """Cramer-Rao lower bound on sd(mu-hat) for the 8-parameter fit."""
    base = _model(p_true, ts)
    Jc = []
    for i in range(8):
        h = 1e-6*max(abs(p_true[i]), 1.0)
        pp = np.array(p_true, float); pp[i] += h
        pm = np.array(p_true, float); pm[i] -= h
        Jc.append(((_model(pp, ts) - _model(pm, ts))/(2*h)).ravel())
    J = np.array(Jc).T
    F = J.T @ J / sigma**2
    C = np.linalg.pinv(F)
    # mu = exp(p6 - p7); d mu = mu (dp6 - dp7)
    g = np.zeros(8); mu = np.exp(p_true[6]-p_true[7]); g[6] = mu; g[7] = -mu
    return float(np.sqrt(g @ C @ g))


SCENARIOS = {
    # name            ball        v0    theta   T(s)   f(Hz)  window
    "pitch":        ("baseball",  40.0,  0.0,   0.46,  300,   30),
    "fly ball":     ("baseball",  45.0, 30.0,   4.60,  100,   40),
    "free kick":    ("soccer",    28.0, 16.0,   1.30,  100,   30),
    "table tennis": ("tabletennis", 12.0, 5.0,  0.30,  500,   40),
}

if __name__ == "__main__":
    print(f"{'scenario':<14}{'L_L (m)':>9}{'path (m)':>10}{'dpsi (deg)':>12}"
          f"{'|w_f|/|w_0|':>13}{'mu true':>9}")
    for name, (ball, v0, thd, T, f, W) in SCENARIOS.items():
        m, r, CD, CL = BALLS[ball]
        kD, kL = coeffs(m, r, CD, CL)
        ts, P, Y = track(v0, np.radians(thd), kD, kL, T, f, z0=2.0)
        s = Y[7, -1]; dpsi = Y[6, -1] - Y[6, 0]
        w0 = np.hypot(Y[3, 0], Y[4, 0]); wf = np.hypot(Y[3, -1], Y[4, -1])
        print(f"{name:<14}{1/kL:>9.2f}{s:>10.2f}{np.degrees(dpsi):>12.3f}"
              f"{wf/w0:>13.4f}{CD/CL:>9.3f}")


def run_e1(name, sigmas, M=4000):
    ball, v0, thd, T, f, W = SCENARIOS[name]
    m, r, CD, CL = BALLS[ball]
    kD, kL = coeffs(m, r, CD, CL)
    ts, P, Y = track(v0, np.radians(thd), kD, kL, T, f, z0=2.0)
    mu = CD/CL
    out = []
    for sg in sigmas:
        est = np.array([E1(ts, P + RNG.normal(0, sg, P.shape), W) for _ in range(M)])
        est = est[np.isfinite(est)]
        out.append((sg, np.mean(est), np.std(est)))
    return mu, len(ts), out
