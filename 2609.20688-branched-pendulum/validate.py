import numpy as np
from sim import run, mode_vars
from params import linear_freqs
# (1) small-amplitude spectra reproduce the linear eigenfrequencies
t, Y = run(1.0, eps=0.005, T=400, n=400001)
dt = t[1]-t[0]
def peaks(x, k):
    X = np.abs(np.fft.rfft((x-x.mean())*np.hanning(len(x)))); f = np.fft.rfftfreq(len(x), dt)
    idx = [i for i in range(1, len(X)-1) if X[i] > X[i-1] and X[i] > X[i+1]]
    idx = sorted(idx, key=lambda i: -X[i])[:k]
    return sorted(round(f[i], 3) for i in idx)
m = mode_vars(Y)
print("theta1 peaks:", peaks(Y[0], 2), " thetaD peak:", peaks(m['tD'], 1), " linear:", ["%.3f" % v for v in linear_freqs()])
# (2) W (their Eq. 13) vs exact dHd/dt, at 47.8 deg, in the early window where theta_D is small
t, Y = run(47.8, T=8, n=80001)
m = mode_vars(Y)
dHd = np.gradient(m['Hd'], t)
for lo, hi in [(0, 2), (2, 3.2), (3.2, 4.2)]:
    s = (t >= lo) & (t < hi)
    r = np.corrcoef(dHd[s], m['W'][s])[0, 1]
    err = np.sqrt(np.mean((dHd[s]-m['W'][s])**2))/np.sqrt(np.mean(dHd[s]**2))
    print(f"t in [{lo},{hi}): max|thD|={np.abs(m['tD'][s]).max():.3f} rad  corr(W,dHd/dt)={r:.4f}  rel rms err={err:.3f}")
