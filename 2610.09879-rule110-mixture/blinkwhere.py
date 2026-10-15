#!/usr/bin/env python3
"""blinkwhere.py -- during a blink, is the extra detection a SPLIT or a new wall?

If the 4th detection sits right beside one of the three real walls, it is one
wall momentarily reported as two, and matching 3 names to the 3 nearest
detections is safe. If it appears in open country, something else is going on.
"""
import numpy as np
from wallcount import step, walls
from tracked import _ring3

N, W, T, EVERY = 300, 15, 4000, 1
rng = np.random.default_rng(777)
import sys
OFFS = {"(1,1,1)": (0,1,2), "(2,2,2)": (0,2,1)}[sys.argv[1] if len(sys.argv)>1 else "(1,1,1)"]
seps, charges, lens = [], [], []
for r in range(6):
    s = _ring3(N, [0, 60, 150], OFFS, rng)
    prev = None
    for t in range(T + 1):
        ch, n = walls(s, t, W)
        p = np.flatnonzero(ch)
        if len(p) > 3 and prev is not None and len(prev) == 3:
            # distance from each extra detection to the nearest previous wall
            d = np.abs(p[:, None] - prev[None, :]); d = np.minimum(d, N - d)
            near = d.min(1)
            extra = np.argsort(near)[-(len(p) - 3):]
            for e in extra:
                seps.append(int(near[e])); charges.append(int(ch[p[e]]))
        if len(p) == 3:
            prev = p
        s = step(s)
seps = np.array(seps); charges = np.array(charges)
print(f"{len(seps)} extra detections over 6 rings x {T} steps "
      f"(= 1 per {6*T/max(1,len(seps)):.0f} steps)")
print(f"distance from the extra detection to the nearest real wall:")
for q in (0, 1, 2, 3, 4, 8, 16, 32, 64):
    print(f"  <= {q:3d} blocks : {(seps <= q).sum():5d}  ({100*(seps<=q).mean():5.1f}%)")
print(f"median {np.median(seps):.0f}, max {seps.max()}")
print(f"charge of the extra detection: "
      f"{dict(zip(*np.unique(charges, return_counts=True)))}")
