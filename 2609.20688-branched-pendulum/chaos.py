"""Is the anti-phase onset the onset of chaos in the in-phase motion? In-phase-only Lyapunov exponent
(two trajectories with thetaD=0 exactly, separated 1e-9 in theta1, renormalised each 1 s) alongside the
anti-phase growth rate (eps=1e-9). Iris."""
import numpy as np, json
from scipy.integrate import solve_ivp
from sim import rhs, run, mode_vars
def lyap(th, T=60.0, d0=1e-9):
    a = np.array([np.radians(th), 0, 0, 0, 0, 0.]); b = a.copy(); b[0] += d0; s = 0.0
    for k in range(int(T)):
        a = solve_ivp(rhs, (0, 1), a, rtol=1e-11, atol=1e-13, method="DOP853").y[:, -1]
        b = solve_ivp(rhs, (0, 1), b, rtol=1e-11, atol=1e-13, method="DOP853").y[:, -1]
        d = np.linalg.norm(b-a); s += np.log(d/d0); b = a + (b-a)*d0/d
    return s/T
def grow(th, T=40.0, eps=1e-9):
    t, Y = run(th, eps=eps, T=T, n=int(T*100)+1); m = mode_vars(Y)
    A = np.log(np.sqrt(m['Hd']/2)); lin = np.abs(m['tD']) < 1e-3
    stop = np.argmax(~lin) if (~lin).any() else len(t)
    return float(np.log(np.abs(m['tD'][:stop]).max()/eps)/t[stop-1]), float(t[stop-1])
def main():
  res = []
  for th in np.arange(38, 50.01, 0.5):
      L = lyap(th); g, tl = grow(th)
      res.append((float(th), L, g, tl)); print(f"{th:5.2f}  in-phase Lyapunov {L:+.3f}/s   anti-phase mean growth {g:+.3f}/s (to {tl:.0f} s)", flush=True)
  json.dump(res, open("chaos.json", "w"))

if __name__ == "__main__":
    main()
