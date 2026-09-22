"""Linear growth rate of the anti-phase mode vs parent angle: eps=1e-9 keeps thetaD linear for tens of
seconds; rate = slope of log(running max |thetaD| amplitude) fitted over the window. Iris."""
import numpy as np, json
from sim import run, mode_vars
eps, T = 1e-9, 40.0
res = []
for th in np.concatenate([np.arange(30, 44, 2.0), np.arange(44, 50.01, 0.25), np.arange(52, 91, 4.0)]):
    t, Y = run(th, eps=eps, T=T, n=int(T*100)+1)
    m = mode_vars(Y)
    A = np.sqrt(m['Hd']/2 + 1e-300)          # amplitude-like, smooth over the oscillation
    lin = np.abs(m['tD']) < 1e-3              # stay in the linear regime
    stop = np.argmax(~lin) if (~lin).any() else len(t)
    tt, la = t[:stop], np.log(A[:stop])
    k = tt > 2
    lam = np.polyfit(tt[k], la[k], 1)[0] if k.sum() > 50 else float('nan')
    res.append((float(th), float(lam), float(tt[-1])))
    print(f"{th:6.2f} deg  lambda = {lam:+.3f} /s   (linear to t={tt[-1]:.1f} s)", flush=True)
json.dump(res, open("rate.json", "w"))
