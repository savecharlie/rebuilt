#!/usr/bin/env python3
"""The detection limit of alpha-hat = 1/CV^2, drawn. Numbers are the measured
output of check.py section 3 (run.log), transcribed, not re-simulated."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

T = [20, 50, 100, 200, 500, 1000]
WRONG = {0.50:[.329,.134,.047,.008,.000,.000], 0.80:[.569,.409,.304,.204,.067,.015],
         0.90:[.625,.525,.438,.359,.235,.141], 1.25:[.195,.207,.172,.116,.043,.009],
         2.00:[.043,.016,.002,.000,.000,.000]}
MED   = {0.50:[.759,.584,.538,.519,.504,.503], 0.80:[1.105,.908,.851,.824,.808,.799],
         0.90:[1.197,1.025,.957,.920,.904,.900], 1.25:[1.582,1.369,1.302,1.272,1.256,1.241],
         2.00:[2.428,2.165,2.088,2.037,1.999,1.985]}
CO = {0.50:"#1b4965", 0.80:"#5fa8d3", 0.90:"#c1666b", 1.25:"#8a9b68", 2.00:"#6b4e71"}

fig, (a, b) = plt.subplots(1, 2, figsize=(11.2, 4.4))
for al, y in WRONG.items():
    a.plot(T, y, "o-", color=CO[al], lw=2, ms=5, label=f"$\\alpha$ = {al}")
a.axhline(.5, color="k", ls=":", lw=1)
a.text(21, .515, "coin flip", fontsize=8.5)
a.set_xscale("log"); a.set_xlabel("record length $T$  (time units, $r$ = 1)")
a.set_ylabel("fraction misclassified across $\\alpha$ = 1")
a.set_title("(a)  how often $\\hat\\alpha=1/CV^2$ lands on the wrong side", fontsize=10.5)
a.legend(fontsize=8.5, frameon=True, framealpha=1.0, edgecolor="none",
         loc="upper right"); a.grid(alpha=.25); a.set_ylim(-0.02, .74)

for al, y in MED.items():
    b.plot(T, y, "o-", color=CO[al], lw=2, ms=5)
    b.axhline(al, color=CO[al], ls="--", lw=.9, alpha=.55)
b.axhline(1.0, color="k", lw=1.4)
b.text(600, 1.03, "SC / QE boundary", fontsize=8.5, ha="right")
b.fill_between([15, 1200], 0.78, 1.0, color="#c1666b", alpha=.07)
b.set_xscale("log"); b.set_xlim(15, 1200); b.set_ylim(0.4, 2.6)
b.set_xlabel("record length $T$  (time units, $r$ = 1)")
b.set_ylabel("median  $\\hat\\alpha$")
b.set_title("(b)  the error is directional: $\\hat\\alpha$ is biased UP", fontsize=10.5)
b.grid(alpha=.25)
b.annotate("true $\\alpha$ = 0.9 reads 1.20 at T=20\nand is still 1.03 at T=50",
           xy=(20, 1.197), xytext=(46, 1.72), fontsize=8.5,
           arrowprops=dict(arrowstyle="->", lw=.9, color="#c1666b"), color="#c1666b")
fig.suptitle("arXiv:2609.30154 — short abundance records systematically call communities STABLE",
             fontsize=11.5, y=.99)
fig.tight_layout(rect=[0, 0, 1, .95])
fig.savefig("detection_limit.png", dpi=155)
print("wrote detection_limit.png")
