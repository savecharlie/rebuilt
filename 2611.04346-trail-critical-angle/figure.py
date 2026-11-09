#!/usr/bin/env python3
"""Two panels for arXiv:2611.04346's equation 2, drawn from the closed forms in verify.py.

(a) the head-wave geometry: ride the trail, leave at theta_c = arccos(C_t/C_g)
(b) the crossover distance d_cross/h against the cost ratio, with its floor at 2

Iris (Opus 5, 1M), 9 Nov 2026.
"""
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

INK, TRAIL, HOT = "#1b1b1b", "#b8860b", "#9b2226"
fig, (ax, bx) = plt.subplots(1, 2, figsize=(12.0, 4.8),
                            gridspec_kw=dict(width_ratios=[1.0, 1.15]))

# ---- (a) geometry, drawn at r = 2 where theta_c is exactly 60 degrees
r = 2.0
h, d = 1.0, 2.6
u = h / math.sqrt(r * r - 1)              # along-trail distance still to run at departure
xq = d - u
ax.axhline(0, color=TRAIL, lw=7, alpha=.30, zorder=0)
ax.axhline(0, color=TRAIL, lw=2.2, zorder=1)
ax.plot([0, xq], [0, 0], color=INK, lw=3.4, zorder=3)
ax.plot([xq, d], [0, h], color=HOT, lw=2.6, zorder=3)
ax.plot([0, d], [0, h], color=INK, lw=1.0, ls=(0, (5, 4)), alpha=.55, zorder=2)
ax.plot([d], [h], marker="o", ms=9, color=INK, zorder=4)
ax.plot([0], [0], marker="o", ms=7, mfc="white", mec=INK, mew=1.6, zorder=4)
ax.plot([xq], [0], marker="o", ms=6, color=HOT, zorder=5)
th = np.linspace(0, math.radians(60), 60)
ax.plot(xq + .30 * np.cos(th), .30 * np.sin(th), color=HOT, lw=1.2)
ax.text(xq + .36, .10, r"$\theta_c$", color=HOT, fontsize=11)
ax.text(d + .55, h * .78, r"$\theta_c=\arccos\!\left(\dfrac{C_{trail}}{C_{ground}}\right)$",
        color=HOT, fontsize=10.5, ha="left", va="center")
ax.annotate("", xy=(d, -.22), xytext=(xq, -.22),
            arrowprops=dict(arrowstyle="<->", color=INK, lw=1.0))
ax.text(d + .10, -.22, r"$u=h/\sqrt{r^{2}-1}$", ha="left", va="center",
        fontsize=10.5, color=INK)
ax.annotate("", xy=(d + .14, 0), xytext=(d + .14, h),
            arrowprops=dict(arrowstyle="<->", color=INK, lw=1.0))
ax.text(d + .24, h / 2, "$h$", va="center", fontsize=11, color=INK)
ax.text(.06, .13, "trail", color=TRAIL, fontsize=11, style="italic")
ax.text(.55, .30, "straight line", color=INK, fontsize=9.5, alpha=.65, rotation=21)
ax.set_xlim(-.30, d + 2.05); ax.set_ylim(-.62, h + .28)
ax.set_aspect("equal"); ax.axis("off")
ax.set_title("a   leave the trail at the critical angle\n"
             r"cost $= C_{trail}\,d + h\sqrt{C_{ground}^2-C_{trail}^2}$"
             "   (drawn at $r=2$, $\\theta_c=60°$)",
             fontsize=11, loc="left", color=INK)

# ---- (b) crossover
rr = np.linspace(1.02, 66, 4000)
dc = 2 * np.sqrt((rr + 1) / (rr - 1))
bx.plot(rr, dc, color=INK, lw=2.2)
bx.axhline(2, color=HOT, lw=1.3, ls=(0, (6, 4)))
bx.text(1.45, 12.6, "floor  $d_{cross}=2h$", color=HOT, fontsize=11)
bx.text(1.45, 11.1, "even a cost-free trail cannot serve two points\n"
                    "closer together than twice their offset",
        color=HOT, fontsize=9.5, va="top")
for rv, lab in ((1.3, "Helbing 30% rule\n$r=1.3$,  $5.54\\,h$"), (2.0, "$r=2$\n$3.46\\,h$")):
    dv = 2 * math.sqrt((rv + 1) / (rv - 1))
    bx.plot([rv], [dv], marker="o", ms=7, color=TRAIL, zorder=5)
    bx.annotate(lab, (rv, dv), textcoords="offset points", xytext=(18, 10),
                fontsize=9.5, color=TRAIL)
bx.set_xscale("log")
bx.set_xlim(1.02, 70); bx.set_ylim(0, 15)
bx.set_xticks([1.25, 2, 3, 5, 10, 20, 40, 65])
bx.set_xticklabels(["1.25", "2", "3", "5", "10", "20", "40", "65"])
bx.minorticks_off()
bx.set_xlabel(r"cost ratio  $r=C_{ground}/C_{trail}$   (Table 1 sweeps 1.25 – 65)", fontsize=10.5)
bx.set_ylabel(r"$d_{cross}/h$", fontsize=11)
bx.set_title("b   how far apart two points must be, along the trail,\n"
             "before the trail is worth using at all", fontsize=11, loc="left", color=INK)
bx.grid(alpha=.18, lw=.7)
for s in ("top", "right"):
    bx.spines[s].set_visible(False)

fig.tight_layout()
fig.savefig("critical_angle_and_crossover.png", dpi=170)
print("wrote critical_angle_and_crossover.png")
