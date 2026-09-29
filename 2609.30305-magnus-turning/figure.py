"""ONE CIRCLE, MANY SPEEDS.

The claim in Gaur 2609.30305 that I found hardest to believe and easiest to draw:
with gravity switched off, a ball spinning about a vertical axis traces a circle
of radius L_L = 2m/(rho C_L A), and DRAG DOES NOT MOVE IT. Not the radius, not the
shape. The drag ball's speed halves thirteen times going round -- the last degree
of the turn takes it most of a day -- and it runs the same ring.

The right panel is the other half of the same algebra. The horizontal velocity,
drawn in the velocity plane, is a logarithmic spiral whose pitch is C_D/C_L and
nothing else: not the mass, not the size, not how hard it was thrown.

Iris, fire 308.
"""
import sys, numpy as np
sys.path.insert(0, ".")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from flight import fly, coeffs, BALLS, turn_event

INK, PAPER, WARM, COOL = "#10100e", "#e8e2d4", "#d4763a", "#4b7d92"

fig = plt.figure(figsize=(13.5, 6.8), facecolor=INK)
gs = fig.add_gridspec(1, 2, wspace=0.14, left=0.055, right=0.975,
                      top=0.855, bottom=0.085)

# ------------------------------------------------------------- left panel
ax = fig.add_subplot(gs[0, 0], facecolor=INK)
m, r, CD, CL = BALLS["baseball"]
kD, kL = coeffs(m, r, CD, CL); LL = 1/kL

sol0 = fly(50.0, 0.0, 0.0, kL, g=0.0, events=[turn_event(2*np.pi)], t_max=200)
T0 = sol0.t_events[0][0]
Y0 = sol0.sol(np.linspace(0, T0, 4000))
ax.plot(Y0[0], Y0[1], color=COOL, lw=9.0, alpha=0.42, solid_capstyle="round",
        zorder=2, label="no drag")

solD = fly(50.0, 0.0, kD, kL, g=0.0, events=[turn_event(2*np.pi)], t_max=6e5)
TD = solD.t_events[0][0]
# sample in TURN ANGLE, not in time -- with drag the last degree takes a day
grid = np.geomspace(1e-3, TD, 60000)
tt = np.interp(np.linspace(0, 2*np.pi, 4000), solD.sol(grid)[6], grid)
YD = solD.sol(tt)
ax.plot(YD[0], YD[1], color=WARM, lw=1.7, zorder=3, solid_capstyle="round",
        label="quadratic drag")

v0 = 50.0; halv = []
for k in range(1, 15):
    v = v0/2**k
    tk = (v0/v - 1)/(kD*v0)          # exact for pure quadratic drag at g = 0
    if tk < TD:
        halv.append((k, solD.sol(tk)))
ax.scatter([h[1][0] for h in halv], [h[1][1] for h in halv], s=26,
           facecolor=INK, edgecolor=WARM, linewidths=1.3, zorder=5)
for k, y in halv[:4]:
    ax.annotate(f"{v0/2**k:.3g} m/s", (y[0], y[1]), textcoords="offset points",
                xytext=(10, -3), color=WARM, fontsize=8.5)
k, y = halv[-1]
ax.annotate(f"{v0/2**k:.3g} m/s", (y[0], y[1]), textcoords="offset points",
            xytext=(16, 20), color=WARM, fontsize=8.5, ha="left",
            arrowprops=dict(arrowstyle="-", color=WARM, alpha=0.5, lw=0.8))

ax.set_aspect("equal")
ax.set_title("both balls run the same ring", color=PAPER, fontsize=13.5,
             loc="left", pad=14)
ax.text(0.50, 0.545, f"radius  $L_L$ = {LL:.2f} m", transform=ax.transAxes,
        color=PAPER, alpha=0.7, fontsize=12, ha="center")
ax.text(0.50, 0.495, "identical to twelve digits", transform=ax.transAxes,
        color=PAPER, alpha=0.45, fontsize=10, ha="center")
ax.text(0.50, 0.415, "rings mark each halving of the", transform=ax.transAxes,
        color=WARM, alpha=0.85, fontsize=9.5, ha="center")
ax.text(0.50, 0.375, "drag ball's speed.  fourteen of them,", transform=ax.transAxes,
        color=WARM, alpha=0.85, fontsize=9.5, ha="center")
ax.text(0.50, 0.335, "and the last degree takes a day", transform=ax.transAxes,
        color=WARM, alpha=0.85, fontsize=9.5, ha="center")
leg = ax.legend(loc="upper right", frameon=False, fontsize=10.5)
for x in leg.get_texts(): x.set_color(PAPER)
print(f"  no drag : one turn in {T0:9.2f} s   50.00 -> "
      f"{np.linalg.norm(sol0.sol(T0)[3:6]):.4f} m/s")
print(f"  drag    : one turn in {TD:9.0f} s   50.00 -> "
      f"{np.linalg.norm(solD.sol(TD)[3:6]):.4f} m/s  "
      f"({TD/86400:.2f} days, {len(halv)} halvings)")

# ------------------------------------------------------------ right panel
ax2 = fig.add_subplot(gs[0, 1], facecolor=INK)
order = sorted(BALLS, key=lambda b: BALLS[b][2]/BALLS[b][3])
cmap = plt.get_cmap("copper")
for i, ball in enumerate(order):
    mm, rr, cd, cl = BALLS[ball]; mu = cd/cl
    psi = np.linspace(0, 2.4*np.pi, 3000)
    w = 50.0*np.exp(-mu*psi)
    col = cmap(0.30 + 0.62*i/(len(order) - 1))
    ax2.plot(w*np.cos(psi), w*np.sin(psi), color=col, lw=2.0,
             solid_capstyle="round")
    ax2.text(0.985, 0.90 - 0.062*i, f"{ball}    $C_D/C_L$ = {mu:.2f}",
             transform=ax2.transAxes, color=col, fontsize=10, ha="right")
ax2.scatter([50.0], [0.0], s=30, color=PAPER, zorder=6, linewidths=0)
ax2.annotate("launch, 50 m/s", (50.0, 0.0), textcoords="offset points",
             xytext=(-12, -17), color=PAPER, alpha=0.6, fontsize=9, ha="right")
ax2.set_aspect("equal"); ax2.set_xlim(-26, 58); ax2.set_ylim(-22, 40)
ax2.set_title("the same flights, drawn in the velocity plane", color=PAPER,
              fontsize=13.5, loc="left", pad=14)
ax2.text(0.015, 0.145, r"$|w| = |w_0|\,e^{-(C_D/C_L)\,\Delta\psi}$",
         transform=ax2.transAxes, color=PAPER, alpha=0.65, fontsize=12.5)
ax2.text(0.015, 0.088, "a logarithmic spiral whose pitch is the drag-to-lift",
         transform=ax2.transAxes, color=PAPER, alpha=0.45, fontsize=9.5)
ax2.text(0.015, 0.046, "ratio, and nothing else",
         transform=ax2.transAxes, color=PAPER, alpha=0.45, fontsize=9.5)

for a in (ax, ax2):
    for sp in a.spines.values(): sp.set_visible(False)
    a.tick_params(colors=PAPER, labelsize=8.5, length=3)
    for lbl in a.get_xticklabels() + a.get_yticklabels(): lbl.set_alpha(0.45)
    a.grid(color=PAPER, alpha=0.07, lw=0.6)

fig.text(0.055, 0.955, "ONE CIRCLE, MANY SPEEDS", color=PAPER, fontsize=17,
         fontweight="bold")
fig.text(0.055, 0.922, "a ball spinning about a vertical axis, with gravity "
         "switched off  ·  after Gaur, arXiv:2609.30305  ·  drawn by Iris",
         color=PAPER, alpha=0.5, fontsize=9.5)
fig.savefig("one_circle_many_speeds.png", dpi=170, facecolor=INK)
print("  wrote one_circle_many_speeds.png")
