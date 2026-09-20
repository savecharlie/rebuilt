import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
r = np.load("rows.npy"); h, p, fx, fc, d, v = r.T
fig, ax = plt.subplots(figsize=(8, 5.2), dpi=130)
ax.semilogy(h, p, "o-", c="#c0392b", label="paper, Table 1 (FFT on a 10 m periodic grid)")
ax.semilogy(h, fc, "x", c="#c0392b", ms=9, label="my rebuild of their FFT, NO trend, NO noise")
ax.semilogy(h, d, "s-", c="#2c3e50", label="same kiln, direct dipole sum (no FFT)")
ax.axhline(4.75, ls="--", c="#7f8c8d")
ax.text(6.6, 3.95, "4.75 nT = mean of the kiln's own field over the 10 m window", color="#555", fontsize=8.5)
lines = ["retention vs the 0.2 m peak,  paper -> direct:"]
for hh in (2, 4, 6, 10):
    i = list(h).index(hh)
    lines.append(f"  {hh:>2} m:  {100*p[i]/p[0]:4.1f}%  ->  {100*d[i]/d[0]:4.2f}%")
ax.text(0.2, 0.3, "\n".join(lines), fontsize=9, family="monospace", va="bottom",
        bbox=dict(fc="white", ec="#bbb"))
ax.set_xlabel("sensor height above ground (m)"); ax.set_ylabel("peak Bz (nT)")
ax.set_title("Hegyi et al. 2609.19349 synthetic kiln: the high-altitude floor is the FFT window's DC term", fontsize=9.5)
ax.legend(fontsize=8.5); ax.grid(alpha=.3, which="both")
fig.tight_layout(); fig.savefig("floor_is_the_window.png")
