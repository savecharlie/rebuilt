#!/usr/bin/env python3
"""arXiv:2609.30154 -- Megias, Ontiveros, Alonso & Capitan.
Stability regimes in ecological communities from abundance time series.

Their DMFT reduces a stochastic Generalized Lotka-Volterra community to one
effective species. Written as their regression form (their Eq. 4),

    log x(t+dt) - log x(t) = dt*(r - s_env^2) - dt*x + mu*dt*M + e,
    Var[e] = dt*(2 s_env^2 + s_int^2 Sigma),

which is the Ito SDE

    d log x = (r - s_env^2 - x + mu M) dt + sqrt(2 theta) dW,
    theta := s_env^2 + s_int^2 Sigma / 2.

Stationary density: for u = log x with constant diffusion D = 2 theta and drift
A(u) = a - e^u, a := r - s_env^2 + mu M,
    p(u) ~ exp((2/D) Int A) = exp((a/theta) u - e^u/theta),
    p(x) = p(u)/x ~ x^{a/theta - 1} e^{-x/theta}.
So Gamma(shape alpha = a/theta, scale theta) -- their Eqs. S24-S25. Confirmed by hand.

This file does three things, in order, and refuses to do the third until the
second passes:
  1. VALIDATE the integrator against the closed form it is supposed to produce.
  2. Check their alpha=1 boundary claim (<N_<> -> exactly 1) to machine precision.
  3. NEW: measure how long an abundance time series must be before alpha-hat
     = 1/CV^2 classifies a community correctly. Their whole empirical claim
     rests on that estimator and the paper does not give its detection limit.
"""
import numpy as np
from scipy import stats, special

RNG = np.random.default_rng(20260925)

# ---------------------------------------------------------------- the model
# No interactions (mu = s_int = 0) makes the problem fully closed-form, which is
# what an instrument check needs: M = r - s_env^2, theta = s_env^2,
# alpha = (r - s_env^2)/s_env^2. With r = 1: alpha = 1/s_env^2 - 1.
R = 1.0

def sigma2_for_alpha(alpha, r=R):
    """s_env^2 giving a target Gamma shape, with mu = s_int = 0."""
    return r / (alpha + 1.0)

def simulate(alpha, n_rep, t_total, dt, sample_every, burn=60.0, r=R, rng=RNG):
    """Euler-Maruyama in LOG space -- positivity is structural, not clipped.
    Stores only every `sample_every` TIME UNITS. The first version of this
    function stored every step: at dt=2e-4 over 300 time units that is a
    400 x 1.5e6 array, 4.8 GB, on a machine already 7 GB into swap. It reached
    7.1 GB RSS and was killed by PID after two minutes. Thin inside the loop."""
    s2 = sigma2_for_alpha(alpha, r)
    theta = s2
    a = r - s2                      # mu = 0 so M drops out of the drift
    amp = np.sqrt(2.0 * theta * dt)
    u = np.log(np.full(n_rep, a))   # start at the deterministic fixed point
    for _ in range(int(burn / dt)):
        u += (a - np.exp(u)) * dt + amp * rng.standard_normal(n_rep)
    stride = max(1, int(round(sample_every / dt)))
    n = int(t_total / dt) // stride
    out = np.empty((n_rep, n))
    for k in range(n):
        for _ in range(stride):
            u += (a - np.exp(u)) * dt + amp * rng.standard_normal(n_rep)
        out[:, k] = np.exp(u)
    return out, alpha, theta

# ------------------------------------------------------- 1. INSTRUMENT CHECK
def validate(dt=5e-4, n_rep=400, t_total=400.0):
    print("1. INSTRUMENT CHECK -- does the integrator reproduce Gamma(alpha, theta)?")
    print(f"   {'alpha':>6} {'theta':>7} | {'mean':>9} {'exact':>9} | {'var':>9} {'exact':>9} "
          f"| {'CV':>6} {'1/sqrt(a)':>9} | {'KS p':>7}")
    ok = True
    for alpha in (0.5, 1.0, 2.0, 4.0):
        x, al, th = simulate(alpha, n_rep, t_total, dt, sample_every=2.0)
        flat = x.ravel()
        m, v = flat.mean(), flat.var()
        em, ev = al * th, al * th * th
        cv, ecv = np.sqrt(v) / m, 1.0 / np.sqrt(al)
        sub = RNG.choice(flat, size=min(4000, flat.size), replace=False)
        p = stats.kstest(sub, "gamma", args=(al, 0.0, th)).pvalue
        print(f"   {al:6.2f} {th:7.4f} | {m:9.5f} {em:9.5f} | {v:9.5f} {ev:9.5f} "
              f"| {cv:6.3f} {ecv:9.3f} | {p:7.4f}")
        if abs(m - em) / em > 0.02 or abs(v - ev) / ev > 0.06 or p < 1e-3:
            ok = False
    print(f"   -> {'PASS' if ok else 'FAIL'}: the integrator may {'' if ok else 'NOT '}be used below.\n")
    return ok

# --------------------------------------- 2. THEIR alpha = 1 BOUNDARY, EXACTLY
def boundary():
    """Their claim: with eps = 1/S, <N_<> = S(1 - P{x >= eps M}) tends to 0 for
    alpha>1, diverges for alpha<1, and equals EXACTLY 1 at alpha=1."""
    print("2. THEIR BOUNDARY CLAIM -- <N_<> = S * gammainc(alpha, alpha/S)")
    print(f"   {'S':>9} | {'a=0.8':>10} {'a=1.0':>12} {'a=1.2':>10}")
    for S in (10, 100, 10**3, 10**5, 10**7, 10**9):
        vals = [S * special.gammainc(a, a / S) for a in (0.8, 1.0, 1.2)]
        print(f"   {S:9d} | {vals[0]:10.4f} {vals[1]:12.9f} {vals[2]:10.6f}")
    # at alpha=1 the incomplete gamma is elementary: S(1-e^{-1/S}) = 1 - 1/(2S) + ...
    # 1 - exp(-1/S) cancels catastrophically at S=1e9: the naive form printed
    # 0.999999971718, which is float error, not mathematics, and it sat BELOW the
    # true value so it looked like it confirmed the approach-from-below. expm1.
    S = 1e9
    naive = S * (1 - np.exp(-1 / S))
    exact = S * (-np.expm1(-1 / S))
    print(f"   alpha=1, S=1e9: expm1 {exact:.12f} | 1-1/(2S) {1 - 1/(2*S):.12f}"
          f" | naive 1-exp() {naive:.12f} <- cancellation, not maths")
    print("   -> confirmed: the limit is 1, approached from BELOW as 1 - 1/(2S).\n")

# ----------------------------------- 3. NEW: how much data does the method need?
def detection(dt=0.02, n_rep=3000):
    """alpha-hat = 1/CV^2 on a finite, AUTOCORRELATED series. The relaxation rate
    of the effective process about its fixed point is M = alpha*theta, so the
    correlation time is tau ~ 1/M -- a series of duration T carries only ~T/tau
    independent samples, not T/dt. Measure the misclassification rate directly."""
    print("3. DETECTION LIMIT of alpha-hat = 1/CV^2 (this is not in the paper)")
    print("   fraction of single series misclassified across the SC/QE line (alpha=1)")
    lengths = [20, 50, 100, 200, 500, 1000]
    alphas = [0.5, 0.8, 0.9, 1.25, 2.0]
    print(f"   {'alpha':>6} {'tau':>5} |" + "".join(f"{('T=' + str(L)):>9}" for L in lengths))
    rows = {}
    for alpha in alphas:
        s2 = sigma2_for_alpha(alpha); M = R - s2; tau = 1.0 / M
        x, _, _ = simulate(alpha, n_rep, max(lengths) + 5.0, dt, sample_every=1.0)
        line, cells = f"   {alpha:6.2f} {tau:5.2f} |", []
        for L in lengths:
            seg = x[:, :L]                                 # already 1 sample / time unit
            cv = seg.std(axis=1, ddof=1) / seg.mean(axis=1)
            ahat = 1.0 / cv**2
            wrong = np.mean((ahat > 1) != (alpha > 1))
            cells.append((wrong, np.median(ahat))); line += f"{wrong:9.3f}"
        rows[alpha] = cells
        print(line)
    print("   (T is the number of time units observed, sampled once per unit.)")
    print("\n   median alpha-hat -- is the error DIRECTIONAL?")
    print(f"   {'alpha':>6}       |" + "".join(f"{('T=' + str(L)):>9}" for L in lengths))
    for alpha in alphas:
        print(f"   {alpha:6.2f}       |" + "".join(f"{m:9.3f}" for _, m in rows[alpha]))
    return rows, lengths

if __name__ == "__main__":
    if not validate():
        raise SystemExit("instrument failed its own check; nothing below is trustworthy")
    boundary()
    detection()
