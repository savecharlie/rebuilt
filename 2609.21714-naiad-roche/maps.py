"""Reproduce their Fig. 1 zero line with their Eq. 1 (MacCullagh), then redraw it with exact
ellipsoid gravity. Shape varied at fixed effective radius (a*b*c fixed at nominal).
Iris, fire 290."""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from roche import surface_g_exact, surface_g_maccullagh, GMp, d, a0, b0, c0

Reff3 = a0 * b0 * c0
def axes(ab, ac):
    a = (Reff3 * ab * ac) ** (1 / 3)
    return a, a / ab, a / ac

def zero_ac(fn, rho, ab):
    f = lambda ac: fn(rho, *axes(ab, ac), GMp, d)
    try: return brentq(f, 1.0001, 4.0)
    except ValueError: return np.nan

for rho in (800., 1280.):
    for ab in (1.0, 1.5, 1.6, 2.0):
        print(f"rho={rho/1000:.2f} a/b={ab}: zero line a/c  MacCullagh {zero_ac(surface_g_maccullagh,rho,ab):.3f}"
              f"   exact {zero_ac(surface_g_exact,rho,ab):.3f}")

ab = np.linspace(1, 2.5, 301); ac = np.linspace(1, 2.5, 301)
AB, AC = np.meshgrid(ab, ac)
fig, axs = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
for ax, rho in zip(axs, (800., 1280.)):
    for fn, col, lab in ((surface_g_maccullagh, "tab:red", "their Eq. 1 (MacCullagh at the surface)"),
                         (surface_g_exact, "tab:blue", "exact ellipsoid gravity")):
        Z = np.vectorize(lambda x, y: fn(rho, *axes(x, y), GMp, d))(AB, AC) * 100
        cs = ax.contour(AB, AC, Z, levels=[0], colors=col, linewidths=2.5)
        ax.plot([], [], color=col, lw=2.5, label=lab)
    ax.errorbar([1.6], [1.846], xerr=[[0.45], [0.45]], yerr=[[0.32], [0.32]], color="k", capsize=4, lw=2)
    ax.set_title(f"rho = {rho/1000:.2f} g/cc: zero net gravity at sub-Neptune point\n(bound below/left of each line)")
    ax.set_xlabel("a/b"); ax.grid(alpha=.3)
axs[0].set_ylabel("a/c"); axs[1].legend(loc="upper right", fontsize=9)
plt.tight_layout(); plt.savefig("naiad_zero_lines.png", dpi=110)
print("wrote naiad_zero_lines.png")
