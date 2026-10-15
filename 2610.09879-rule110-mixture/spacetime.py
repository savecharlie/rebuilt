#!/usr/bin/env python3
"""spacetime.py -- LOOK at the wall before believing any velocity.

veltri.py reports wall velocities that change sign and go to zero across the 32
implementations.  Those numbers come from a detector built and validated for
(36,37,5); a wall that is wide, chaotic, or absent would also produce them, and
the detector cannot tell me which.  So: draw the raw phase field in space-time
for a handful of triples and read the picture.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

from tri import Kit

N, T, W = 600, 700, 15
PICK = [(36, 37, 5), (38, 37, 21), (52, 37, 5), (52, 37, 29),
        (54, 37, 13), (36, 39, 21), (54, 37, 5), (38, 39, 5)]

cmap = ListedColormap(["#1b1b22", "#d8552a", "#2f8fd0", "#e8d44d"])  # -1,0,1,2

fig, axes = plt.subplots(2, 4, figsize=(19, 10), sharex=True, sharey=True)
for ax, tri in zip(axes.ravel(), PICK):
    k = Kit(tri)
    s = k.two_domain(N, 1, np.random.default_rng(917))
    img = np.zeros((T, N), np.int8)
    amb = np.zeros(T)
    for t in range(T):
        lab = k.labels(s, t)
        img[t] = lab
        amb[t] = (lab < 0).mean()
        s = k.step(s)
    ax.imshow(img + 1, cmap=cmap, aspect="auto", interpolation="nearest",
              vmin=0, vmax=3, origin="upper")
    ax.set_title(f"{tri}   unlabelled {amb.mean()*100:.0f}%", fontsize=11)
    ax.set_xlabel("block")
axes[0, 0].set_ylabel("step")
axes[1, 0].set_ylabel("step")
fig.suptitle("RAW beat-phase label, two-domain ring, same initial condition.  "
             "dark = unlabelled (form `100`, legal in all three beats); "
             "orange/blue/yellow = beats 0/1/2", fontsize=13)
fig.tight_layout()
fig.savefig("spacetime_by_triple.png", dpi=110)
print("wrote spacetime_by_triple.png")

# and the quantitative companion: how much of the ring is NOT a clean domain
print("\nfraction of blocks that are unlabelled, far from the two built walls")
print("(a locked domain is ~41% unlabelled -- that is the form `100`; much more "
      "means the domain itself has stopped being locked)")
for tri in PICK:
    k = Kit(tri)
    s = k.two_domain(N, 1, np.random.default_rng(917))
    for _ in range(300):
        s = k.step(s)
    lab = k.labels(s, 300)
    # the two built walls started at 150 and 450; look in the middle of each domain
    mid = np.r_[lab[280:320], lab[20:60]]
    print(f"   {str(tri):<16} whole ring {np.mean(lab<0)*100:5.1f}%   "
          f"domain interiors {np.mean(mid<0)*100:5.1f}%")
