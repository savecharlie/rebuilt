"""LAW 1: validate the instrument against answers I already know, BEFORE it is
allowed to say anything new.

Three closed forms that have nothing to do with the paper:
  A. vacuum parabola          (kD = kL = 0)
  B. 1-D fall with quadratic drag, v = v_t tanh(g t / v_t)   (kL = 0, vertical)
  C. horizontal coast under pure quadratic drag, no gravity, no lift

If any of these is not at the integrator's tolerance, nothing below it means
anything.
"""
import numpy as np
from flight import fly, rhs, G, ground_event
from scipy.integrate import solve_ivp

fails = []
def check(name, got, want, tol):
    rel = abs(got - want) / max(abs(want), 1e-300)
    ok = rel < tol
    print(f"  [{'ok ' if ok else 'FAIL'}] {name:<46s} got {got:.15g}  want {want:.15g}   rel {rel:.2e}")
    if not ok: fails.append(name)

print("A. vacuum parabola (kD = kL = 0)")
for v0, thdeg in [(50, 30), (30, 60), (120, 15), (10, 80)]:
    th = np.radians(thdeg)
    sol = fly(v0, th, 0.0, 0.0, t_max=100, events=[ground_event()])
    t = sol.t_events[0][0]; y = sol.y_events[0][0]
    check(f"flight time v0={v0} th={thdeg}", t, 2*v0*np.sin(th)/G, 1e-11)
    check(f"range      v0={v0} th={thdeg}", np.hypot(y[0], y[1]),
          v0*v0*np.sin(2*th)/G, 1e-11)
    check(f"path length v0={v0} th={thdeg}", y[7],
          (v0*v0/G)*(np.sin(th) + np.cos(th)**2*np.log(np.tan(th)+1/np.cos(th))), 1e-11)

print("\nB. 1-D fall, quadratic drag (analytic tanh)")
kD = 0.01
vt = np.sqrt(G/kD)
for T in [1.0, 5.0, 20.0]:
    # straight down-start from rest: integrate the vertical component alone
    f = lambda t, y: [y[1], -G - kD*abs(y[1])*y[1]]
    s = solve_ivp(f, (0, T), [0.0, 0.0], method="DOP853", rtol=1e-13, atol=1e-14)
    z, v = s.y[0][-1], s.y[1][-1]
    check(f"v(t={T})", v, -vt*np.tanh(G*T/vt), 1e-10)
    check(f"z(t={T})", z, -(vt*vt/G)*np.log(np.cosh(G*T/vt)), 1e-10)

print("\nC. horizontal coast, quadratic drag, g = 0, kL = 0")
kD = 0.004; v0 = 50.0
for T in [0.5, 4.0, 30.0]:
    sol = solve_ivp(rhs, (0, T), [0,0,0, v0,0,0, 0,0], args=(kD, 0.0, +1, 0.0),
                    method="DOP853", rtol=1e-13, atol=1e-14)
    y = sol.y[:, -1]
    check(f"v(t={T})", y[3], v0/(1 + kD*v0*T), 1e-11)
    check(f"x(t={T})", y[0], np.log1p(kD*v0*T)/kD, 1e-11)
    check(f"s(t={T})", y[7], np.log1p(kD*v0*T)/kD, 1e-11)

print()
if fails:
    print(f"INSTRUMENT NOT VALIDATED -- {len(fails)} failures:", fails)
    raise SystemExit(1)
print("instrument validated on 3 independent closed forms. Proceed.")
