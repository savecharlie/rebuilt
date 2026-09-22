import numpy as np, json, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from sim import run, mode_vars
from wclosed import W_closed
from chaos import grow
fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
t, Y = run(47.8, eps=1e-4, T=1.5, n=30001); m = mode_vars(Y); dHd = np.gradient(m['Hd'], t); Wc = W_closed(m, Y)[0]
ax[0].plot(t, dHd*1e8, c="k", lw=2.4, label="exact dH_d/dt")
ax[0].plot(t, m['W']*1e8, c="#c43", lw=1.2, label="paper Eq. 13 (W_K + W_P)")
ax[0].plot(t, Wc*1e8, c="#3a7", lw=1, ls="--", label="+ W_X + W_I (this rebuild)")
ax[0].set_xlabel("t − t0 (s)"); ax[0].set_ylabel("×1e-8 J/s"); ax[0].legend(fontsize=8)
ax[0].set_title("θ1(0) = 47.8°, θ_D(0) = 1e-4")
eps = np.array([1e-2, 1e-3, 1e-4, 1e-5]); ep, ec = [], []
for e in eps:
    t, Y = run(47.8, eps=e, T=1.5, n=30001); m = mode_vars(Y); d = np.gradient(m['Hd'], t); s = slice(200, -200)
    r = lambda W: np.sqrt(np.mean((d[s]-W[s])**2))/np.sqrt(np.mean(d[s]**2))
    ep.append(r(m['W'])); ec.append(r(W_closed(m, Y)[0]))
ax[1].loglog(eps, ep, "o-", c="#c43", label="paper Eq. 13"); ax[1].loglog(eps, ec, "s--", c="#3a7", label="corrected")
ax[1].loglog(eps, 1.4*(eps/1e-2)**2*1e-4/1.4, ":", c="gray", label="∝ ε²")
ax[1].set_xlabel("initial disagreement ε (rad)"); ax[1].set_ylabel("rel. rms error vs exact"); ax[1].legend(fontsize=8)
ax[1].set_title("a leading-order formula must converge")
L = {}
for f in __import__("glob").glob("lyap_*.log"):
    a = open(f).read().split()
    if len(a) >= 3: L[float(a[0])] = float(a[2].split(":")[1])
ths = sorted(k for k in L if 38 <= k <= 50)
ax[2].plot(ths, [L[k] for k in ths], "o-", c="#335", label="in-phase Lyapunov (240 s)")
gt = [41, 42, 42.25, 42.5, 42.75, 43, 43.5, 44, 46]
ax[2].plot(gt, [grow(x)[0] for x in gt], "s-", c="#c43", label="anti-phase growth rate")
ax[2].axvspan(47.6, 47.8, color="gold", alpha=.5, label="experiment 47.6–47.8°")
ax[2].set_xlabel("θ1(0) (deg)"); ax[2].set_ylabel("/s"); ax[2].legend(fontsize=8); ax[2].set_title("both switch on at 42.25–42.5°")
plt.tight_layout(); plt.savefig("branched_pendulum_rebuild.png", dpi=85)
