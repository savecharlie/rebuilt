import numpy as np
from sim import run, mode_vars
for th in (20.0, 47.8):
  for eps in (1e-2, 1e-3, 1e-4, 1e-5):
    t, Y = run(th, eps=eps, T=1.5, n=30001)
    m = mode_vars(Y)
    dHd = np.gradient(m['Hd'], t)
    s = slice(200, -200)
    err = np.sqrt(np.mean((dHd[s]-m['W'][s])**2))/np.sqrt(np.mean(dHd[s]**2))
    fit = np.dot(m['W'][s], dHd[s])/np.dot(m['W'][s], m['W'][s])
    print(f"theta1(0)={th:5.1f}  eps={eps:.0e}  max|thD|={np.abs(m['tD']).max():.1e}  rel rms err of W = {err:.3f}   best scale dHd/W = {fit:.3f}")
