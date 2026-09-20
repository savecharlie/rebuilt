import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from parse import read_ndk

e = read_ndk("jan76_dec20.ndk"); M0 = e["M0"]*1e-7
thr = e["Mrr"] > 0
d, m = e["dep"][thr], M0[thr]
o = np.argsort(d); d, m = d[o], m[o]
cum = 100*np.cumsum(m)/m.sum()

A_dot, h, drho, g, YR = 3.0e6, 100e3, 50.0, 9.8, 365.25*86400
zz = np.linspace(0, 2890, 600)
P  = drho*g*(zz*1e3)*A_dot*h/YR/1e12

fig, ax = plt.subplots(1, 2, figsize=(11.5, 6.0))

ax[0].plot(cum, d, lw=2.2, color="#1b3a5c")
ax[0].axhspan(0, 70, color="#c0392b", alpha=.10)
ax[0].axhline(70, color="#c0392b", ls="--", lw=1.2)
ax[0].text(4, 150, "95.0% of the thrust moment\nis shallower than 70 km\n(moment-weighted mean 38 km)",
           color="#c0392b", fontsize=9.5)
ax[0].set_xlabel("cumulative % of thrust ($M_{rr}>0$) seismic moment")
ax[0].set_ylabel("centroid depth (km)")
ax[0].set_title("Where the thrust earthquakes are\nGCMT 1976–2020, 29,846 events", fontsize=10.5)
ax[0].set_ylim(700, 0); ax[0].set_xlim(0, 100); ax[0].grid(alpha=.25)

ax[1].plot(P, zz, lw=2.2, color="#1b3a5c")
ax[1].axhspan(0, 70, color="#c0392b", alpha=.10)
ax[1].axhline(70, color="#c0392b", ls="--", lw=1.2)
ax[1].axhline(660, color="#888", ls=":", lw=1.0)
ax[1].text(14.6, 645, "660 km", fontsize=8, color="#666", ha="right")
ax[1].axvspan(10, 14, color="#e8a33d", alpha=.30)
ax[1].axvspan(0.3, 1.3, color="#4a8f4a", alpha=.45)
ax[1].plot([0.33], [70], "o", color="#c0392b", zorder=5)
ax[1].annotate("0.33 TW has been released\nby the time the slab is 70 km down",
               (0.33, 70), xytext=(2.3, 420), fontsize=9.5, color="#c0392b",
               arrowprops=dict(arrowstyle="->", color="#c0392b", lw=1.0))
ax[1].text(12.0, 1450, "10–14 TW\n\nthe figure the paper\nassigns to these\nearthquakes — and the\ndepth you must reach\nto have released it",
           fontsize=9.5, color="#7a4f0c", ha="center", va="center")
ax[1].text(0.8, 2450, "0.3–1.3 TW\nwhat island-arc magma\nproduction requires",
           fontsize=9.5, color="#255225", ha="left", va="center")
ax[1].set_xlabel("gravitational power released by slab descent through the top $z$ km  (TW)")
ax[1].set_ylabel("depth (km)")
ax[1].set_title("Where the gravitational energy is\n3 km²/yr of plate × 100 km thick × Δρ 50 kg/m³", fontsize=10.5)
ax[1].set_xlim(0, 15); ax[1].set_ylim(2890, 0); ax[1].grid(alpha=.25)

fig.suptitle("arXiv:2609.18036 — the conversion is placed where the energy is not", fontsize=12.5)
fig.tight_layout()
fig.savefig("where_the_moment_is.png", dpi=150)
print("ok")
