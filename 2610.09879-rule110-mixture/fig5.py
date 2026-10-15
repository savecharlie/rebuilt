#!/usr/bin/env python3
"""fig5.py -- the figure for section 5: one rule 110, 32 implementations,
three kinds of wall, and a self-organisation result that does not survive.
"""
import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

from tri import MAIN_CUBE, Kit
from bits import word

SW = json.load(open("sweep32b.json"))["data"]
LKN = 128
try:
    LK = json.load(open("lock32_N128.json"))["data"]
except FileNotFoundError:
    # lock32.py writes its JSON only at the end and takes hours on the
    # implementations whose rings never lock (every run burns the full budget).
    # Its per-triple lines are already on stdout, so read those.
    import re
    LK = {}
    for ln in open("lock32.out"):
        m = re.match(r"\s+(\(\d+, \d+, \d+\))\s+locked ([\d.]+) \+- ([\d.]+)", ln)
        if m:
            LK[m.group(1)] = {"frac": float(m.group(2)), "se": float(m.group(3))}
    print(f"   (read {len(LK)} locking fractions from lock32.out)")
R = json.load(open("sweep32b.json"))["R"]

SHOW = [(36, 37, 5), (36, 37, 13), (52, 37, 29), (38, 37, 21)]
LABEL = ["sharp mover\n(the paper's triple)",
         "INVADER: a phase legal\non no beat, growing",
         "third-beat sliver,\npinned", "sharp mover, faster"]
# five categories: the label field's -1 has TWO causes and they are opposite
cmap = ListedColormap(["#d12f4f",   # `001` -- legal on NO beat
                       "#33333d",   # `100` -- legal on every beat
                       "#e0603a", "#3f8fc9", "#f2e34a"])   # beats 0/1/2

fig = plt.figure(figsize=(17, 12))
gs = fig.add_gridspec(3, 4, height_ratios=[1.35, 0.85, 0.9], hspace=0.55, wspace=0.26)

# --- row 1: space-time ------------------------------------------------
N, T = 400, 460
for j, (tri, lab) in enumerate(zip(SHOW, LABEL)):
    ax = fig.add_subplot(gs[0, j])
    k = Kit(tri)
    s = k.two_domain(N, 1, np.random.default_rng(31))
    img = np.zeros((T, N), np.int8)
    for t in range(T):
        lfield = k.labels(s, t)          # NOT `lab` -- that is the panel label
        b = s.reshape(-1, 3).astype(np.int64)
        code = 4 * b[:, 0] + 2 * b[:, 1] + b[:, 2]
        img[t] = np.where(lfield >= 0, lfield + 2, np.where(code == 1, 0, 1))
        s = k.step(s)
    ax.imshow(img, cmap=cmap, aspect="auto", interpolation="nearest", vmin=0, vmax=4)
    f = LK.get(str(tri), {}).get("frac")
    ft = f"locks {f:.0%} of random rings" if f is not None else "locking not yet measured"
    ax.set_title(f"{tri}\n{lab}\n{ft}",
                 fontsize=10.5, linespacing=1.4)
    ax.set_xticks([]); ax.set_yticks([0, 200, 400])
    if j == 0:
        ax.set_ylabel("step")

# --- row 2 left: the core profiles ------------------------------------
ax = fig.add_subplot(gs[1, :2])
off = np.arange(-R, R + 1)
for tri, lab, col in zip(SHOW, LABEL, ["#2b6cb0", "#c05621", "#b7791f", "#2f855a"]):
    e = SW[str(tri)]
    if e.get("v1") is None:
        continue
    ax.plot(off, e["mid1"], lw=2, color=col, label=f"{tri}  MID")
    ax.plot(off, np.array(e["amb1"]) - np.r_[e["amb1"][:10], e["amb1"][-10:]].mean(),
            lw=1.1, ls="--", color=col, alpha=0.8)
ax.axhline(0, color="k", lw=0.5)
ax.set_xlabel("blocks from the wall centre")
ax.set_ylabel("occupancy")
ax.set_title("inside a charge-1 wall: solid = third beat, dashed = excess "
             "unlabelled blocks over bulk\n(the dashed channel mixes the all-beats "
             "form with the no-beat one -- see invade.py)", fontsize=9.5)
ax.legend(fontsize=8.5, ncol=2)
ax.set_xlim(-30, 30)

# --- row 2 right: core vs velocity ------------------------------------
ax = fig.add_subplot(gs[1, 2:])
for c, mk, col in ((1, "o", "#2b6cb0"), (2, "s", "#c05621")):
    x, y = [], []
    for t in MAIN_CUBE:
        e = SW[str(t)]
        if e.get(f"v{c}") is None:
            continue
        x.append(max(e[f"mid{c}"])); y.append(abs(e[f"v{c}"][0]))
    ax.scatter(x, y, marker=mk, s=46, color=col, alpha=0.75, label=f"charge {c}")
ax.axhline(0.01, color="k", lw=0.6, ls=":")
ax.text(0.45, 0.013, "pinned below here", fontsize=8.5)
ax.set_xlabel("third-beat peak occupancy inside the wall")
ax.set_ylabel("|wall velocity|  (blocks/step)")
ax.set_title("the core does not predict the speed with one sign:\n"
             "charge-1 walls pin WITH a big core, charge-2 walls pin WITHOUT one",
             fontsize=10)
ax.legend(fontsize=9)

# --- row 3: locking fraction by implementation ------------------------
ax = fig.add_subplot(gs[2, :])
HAVE = [t for t in MAIN_CUBE if str(t) in LK]
order = sorted(HAVE, key=lambda t: -LK[str(t)]["frac"])
f = [LK[str(t)]["frac"] for t in order]
se = [LK[str(t)]["se"] for t in order]
cols = ["#2f855a" if LK[str(t)]["frac"] > 0.5 else "#b03030" for t in order]
ax.bar(range(len(order)), f, yerr=se, color=cols, edgecolor="k", lw=0.4, capsize=2)
ax.set_xticks(range(len(order)))
ax.set_xticklabels([f"{t[0]},{t[1]},{t[2]}" for t in order], rotation=90, fontsize=7.5)
ax.axhline(0.942, color="k", ls="--", lw=1)
ax.text(len(order)-0.6, 0.885, "the paper's triple: 0.942, rank 1 of 32  (median of the 32 is 0.083)", ha="right", fontsize=9)
ax.set_ylabel("fraction of random rings that\nlock within 60,000 steps  (N = 128)")
ax.set_title("all 32 emulate rule 110 exactly, step for step, with the same encoding "
             "— they disagree only on five table entries no computation ever reaches",
             fontsize=11)
ax.set_ylim(0, 1.0)

fig.suptitle("One rule 110, thirty-two implementations: the defect is where the choice shows\ncrimson = the block form legal on NO beat, charcoal = the form legal on every beat",
             fontsize=14, y=0.99)
fig.savefig("section5_implementations.png", dpi=108, bbox_inches="tight")
print("wrote section5_implementations.png")
