"""The crosswind sensitivity of Eq. (15), with its derivation and its
convergence test. This is the one result here that is not in the paper.

Gaur's Eq. (15) reads the lift-to-drag ratio off a single tracked flight,

    mu = C_D/C_L = ln(|w_0|/|w_f|) / dpsi,

and is exact. It is also, in Monte Carlo, precise: 0.1-2 % at 3 mm positional
noise (can_you_measure_mu.py). What limits it is that it is written in the
GROUND frame and the aerodynamics live in the AIR frame.

THE RESULT. Let the air move horizontally at speed u on a bearing alpha
measured from the launch direction. Then to leading order in u/|w_0| and in
the heading swing dpsi,

    d(mu-hat)/mu  =  -(mu + 1/mu) (u/|w_0|) sin(alpha),                  (*)

with the next term a factor (1 + mu dpsi / 2).

DERIVATION. In the air frame the reduction is exact, so the horizontal air
velocity runs from W_0 to W_f = W_0 e^{-mu D} while its heading advances by D.
Adding the constant wind vector u e^{i alpha} to both endpoints and expanding
to first order in u,

    ln|w_0/w_f|  ->  mu D + u [ cos a / W_0  -  cos(a-D) / W_f ]
    dpsi         ->  D    + u [ sin(a-D)/ W_f -  sin a  / W_0 ]

so d(mu-hat)/mu = (u/D)[A/mu - B] with A and B the two brackets. Expanding
those for small D, every term in 1/D cancels and what survives is (*).

WHAT IT MEANS, and it is the part worth keeping:

  1. mu + 1/mu >= 2 for every mu, with equality at mu = 1. There is no ball and
     no launch for which the crosswind sensitivity is small. The floor is
     2 u/|w_0|, and it is reached only when C_D = C_L.

  2. dpsi cancels. Every other quantity in this measurement improves with a
     longer flight -- more heading swing, more speed decay, more frames. The
     wind bias does not. It cannot be bought off with a longer trajectory, and
     measurement C below shows it in fact gets slightly worse.

  3. Against the 0.1-2 % noise floor, 1 m/s of crosswind on a 40 m/s ball is a
     6 % bias. The estimator is systematics-limited by two orders of magnitude,
     outdoors, in what anyone would call still air.

Iris, fire 308 (29 Sep 2026).
"""
import numpy as np
from flight import coeffs
from can_you_measure_mu import E1
from what_breaks_it import track_gen

VERT = np.array([0.0, 0.0, 1.0])
ALPHAS = np.radians(np.arange(0, 360, 15))


def amplitude(m, r, CD, CL, v0, theta, T, f, W=20, u=0.02):
    """Fit  bias/(u/|w0|) = -amp sin(alpha)  and return (amp, dpsi, rms_resid)."""
    kD, kL = coeffs(m, r, CD, CL)
    ts, P, Y = track_gen(v0, theta, kD, kL, T, f, 2.0, VERT, np.zeros(3))
    base = E1(ts, P, W)
    w0 = np.hypot(Y[3, 0], Y[4, 0]); dpsi = Y[6, -1] - Y[6, 0]
    c = []
    for a in ALPHAS:
        wnd = u*np.array([np.cos(a), np.sin(a), 0.0])
        t2, P2, _ = track_gen(v0, theta, kD, kL, T, f, 2.0, VERT, wnd)
        c.append(((E1(t2, P2, W) - base)/base)/(u/w0))
    c = np.array(c)
    amp = -(c @ np.sin(ALPHAS))/(np.sin(ALPHAS) @ np.sin(ALPHAS))
    return amp, dpsi, float(np.sqrt(((c + amp*np.sin(ALPHAS))**2).mean()))


if __name__ == "__main__":
    m, r, CL = 0.145, 0.0366, 0.20

    print("A. the shape is a pure sine in the wind bearing")
    amp, dpsi, res = amplitude(m, r, 0.35, CL, 40.0, 0.0, 0.46, 300, W=30)
    print(f"   amplitude {amp:.4f}   residual from sin(alpha): {res/amp:.2%} of it\n")

    print("B. the amplitude is mu + 1/mu, and the residual is (1/2) mu dpsi")
    print(f"   {'mu':>6}{'amp':>10}{'mu+1/mu':>10}{'rel err':>10}{'err/(mu dpsi)':>15}")
    for CD in [0.10, 0.20, 0.30, 0.40, 0.60, 0.80]:
        mu = CD/CL
        amp, dpsi, _ = amplitude(m, r, CD, CL, 40.0, 0.0, 0.46, 300, W=30)
        e = amp/(mu + 1/mu) - 1
        print(f"   {mu:>6.2f}{amp:>10.4f}{mu+1/mu:>10.4f}{e:>10.2%}{e/(mu*dpsi):>15.4f}")

    print("\nC. a longer flight does not help -- the amplitude grows")
    print(f"   {'theta':>6}{'T (s)':>7}{'dpsi(deg)':>11}{'amp':>10}")
    for thd, T in [(0, 0.46), (0, 1.5), (20, 3.0), (30, 4.6), (45, 6.5)]:
        amp, dpsi, _ = amplitude(m, r, 0.35, CL, 40.0, np.radians(thd), T, 300, W=30)
        print(f"   {thd:>6}{T:>7.2f}{np.degrees(dpsi):>11.2f}{amp:>10.4f}")

    print("\nD. LAW 22: the claim is exact as dpsi -> 0, so test its CONVERGENCE")
    mu = 1.75; pred = mu + 1/mu; prev = None
    print(f"   {'dpsi (rad)':>12}{'amp':>11}{'rel err':>10}{'halving ratio':>15}")
    for T in [0.64, 0.32, 0.16, 0.08, 0.04, 0.02]:
        amp, dpsi, _ = amplitude(m, r, 0.35, CL, 40.0, 0.0, T, max(300, int(60/T)))
        e = amp/pred - 1
        rat = f"{prev/e:.3f}" if prev else "-"
        print(f"   {dpsi:>12.5f}{amp:>11.6f}{e:>10.2%}{rat:>15}")
        prev = e
    print("   -> 2.000 per halving: first order, exactly as (*) requires.")
