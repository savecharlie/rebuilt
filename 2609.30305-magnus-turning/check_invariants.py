"""The four turning invariants of Gaur 2609.30305, checked against my own
integration of Eq. (1). Nothing here uses the paper's reduced equation.
"""
import numpy as np
from flight import fly, coeffs, G, BALLS, turn_event, ground_event

def L(kL): return 1.0/kL

rows = []
def rec(name, dev, note=""):
    rows.append((name, dev, note))
    print(f"  {name:<52s} max dev {dev:.2e}   {note}")

print("=== IV.1  dpsi/ds = kL, pointwise along the path ===")
for ball in ["baseball", "frisbee", "tabletennis"]:
    m, r, CD, CL = BALLS[ball]
    for drag in (False, True):
        kD, kL = coeffs(m, r, CD if drag else 0.0, CL)
        worst = 0.0
        for v0, thd in [(50, 30), (30, 60), (120, 15), (15, 75)]:
            sol = fly(v0, np.radians(thd), kD, kL, events=[ground_event()], t_max=200)
            ts = np.linspace(0, sol.t_events[0][0], 400)[1:-1]
            Y = sol.sol(ts)
            # dpsi/ds = psidot / |v|, both recomputed from the state
            V = np.linalg.norm(Y[3:6], axis=0)
            # psidot from the stored psi channel, by differentiating the ODE again
            vx, vy, vz = Y[3], Y[4], Y[5]
            ax = -kD*V*vx + kL*V*(-vy); ay = -kD*V*vy + kL*V*(vx)
            psidot = (vx*ay - vy*ax)/(vx*vx + vy*vy)
            worst = max(worst, np.max(np.abs(psidot/(kL*V) - 1.0)))
        rec(f"{ball}, drag {'on ' if drag else 'off'}", worst)

print("\n=== IV.1  a 180 deg turn costs exactly pi*L_L ===")
# Use the frisbee, which the paper says can actually close a half turn.
m, r, CD, CL = BALLS["frisbee"]
kD, kL = coeffs(m, r, CD, CL)
target = np.pi*L(kL)
worst = 0.0; n = 0
for v0 in np.linspace(10, 120, 12):
    for thd in [5, 20, 45, 60, 80]:
        sol = fly(v0, np.radians(thd), kD, kL, events=[turn_event(np.pi)], t_max=600)
        if len(sol.t_events[0]) == 0:   # never made the half turn
            continue
        s = sol.y_events[0][0][7]; n += 1
        worst = max(worst, abs(s/target - 1))
rec(f"S(180 deg) / (pi L_L) - 1, drag ON, {n} launches", worst,
    f"pi L_L = {target:.6f} m")

print("\n=== IV.2  |w| = |w0| exp(-mu * dpsi) ===")
for ball in BALLS:
    m, r, CD, CL = BALLS[ball]
    kD, kL = coeffs(m, r, CD, CL); mu = CD/CL
    sol = fly(60.0, np.radians(40), kD, kL, events=[ground_event()], t_max=300)
    ts = np.linspace(0, sol.t_events[0][0], 300)[1:]
    Y = sol.sol(ts)
    w = np.hypot(Y[3], Y[4]); w0 = np.hypot(*fly(60.0, np.radians(40), kD, kL,
                                                dense=False, t_max=1e-9).y[3:5, 0])
    dev = np.max(np.abs(w/(w0*np.exp(-mu*Y[6])) - 1))
    rec(f"{ball} (mu = {mu:.3f})", dev)

print("\n=== IV.3  R_h = L_L cos(gamma) ===")
for ball in ["baseball", "soccer", "frisbee"]:
    m, r, CD, CL = BALLS[ball]
    for drag in (False, True):
        kD, kL = coeffs(m, r, CD if drag else 0.0, CL)
        sol = fly(55.0, np.radians(50), kD, kL, events=[ground_event()], t_max=300)
        ts = np.linspace(0, sol.t_events[0][0], 400)[1:-1]
        Y = sol.sol(ts)
        vx, vy, vz = Y[3], Y[4], Y[5]
        V = np.sqrt(vx*vx+vy*vy+vz*vz); w = np.hypot(vx, vy)
        ax = -kD*V*vx + kL*V*(-vy); ay = -kD*V*vy + kL*V*vx
        psidot = (vx*ay - vy*ax)/(w*w)
        Rh = w/psidot                      # (ds_h/dt)/(dpsi/dt)
        gamma = np.arctan2(vz, w)
        rec(f"{ball}, drag {'on ' if drag else 'off'}",
            np.max(np.abs(Rh/(L(kL)*np.cos(gamma)) - 1)))

print("\n=== IV.4  g = 0 : the ground track is a circle of radius L_L ===")
def fit_circle(x, y):
    A = np.c_[2*x, 2*y, np.ones_like(x)]
    b = x*x + y*y
    c, *_ = np.linalg.lstsq(A, b, rcond=None)
    return np.sqrt(c[2] + c[0]**2 + c[1]**2)
for ball in ["baseball", "frisbee"]:
    m, r, CD, CL = BALLS[ball]
    for drag in (False, True):
        kD, kL = coeffs(m, r, CD if drag else 0.0, CL)
        sol = fly(50.0, 0.0, kD, kL, g=0.0, events=[turn_event(np.pi)], t_max=4000)
        T = sol.t_events[0][0]
        ts = np.linspace(0, T, 2000)
        Y = sol.sol(ts)
        R = fit_circle(Y[0], Y[1])
        vend = np.linalg.norm(sol.y_events[0][0][3:6])
        rec(f"{ball}, drag {'on ' if drag else 'off'}", abs(R/L(kL) - 1),
            f"R={R:.9f}  L_L={L(kL):.9f}  speed 50 -> {vend:.3f}")

print("\n=== summary ===")
worst = max(r[1] for r in rows)
print(f"  {len(rows)} checks, largest relative deviation {worst:.2e}")
