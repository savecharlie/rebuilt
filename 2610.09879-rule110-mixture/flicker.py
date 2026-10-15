#!/usr/bin/env python3
"""flicker.py -- is a wall-count change a physical event or a one-frame blink?

measure4b.py's hard stop ended (1,1,1) tracks after a median of 650 steps, while
fire 341 reported a median resolution of 54,900 steps for the same rings. A
factor of 80 between two of my own measurements means at least one is about the
detector. So: do not stop. Record the raw count every frame and LOOK at it.
"""
import numpy as np
from wallcount import step, walls, locked_ring
from tracked import _ring3

N, W, T, EVERY = 300, 15, 6000, 5

for name, offs in (("(1,1,1)", (0, 1, 2)), ("(2,2,2)", (0, 2, 1))):
    rng = np.random.default_rng(777)
    print(f"\n=== {name} ===")
    for r in range(4):
        s = _ring3(N, [0, 60, 150], offs, rng)
        cnt, ts = [], []
        for t in range(T + 1):
            if t % EVERY == 0:
                _, n = walls(s, t, W)
                cnt.append(int(n)); ts.append(t)
            s = step(s)
        c = np.array(cnt)
        first = np.flatnonzero(c != c[0])
        runs = []
        i = 0
        while i < len(c):
            j = i
            while j + 1 < len(c) and c[j + 1] == c[i]:
                j += 1
            runs.append((c[i], ts[i], ts[j]))
            i = j + 1
        comp = " ".join(f"{v}@{a}-{b}" for v, a, b in runs[:14])
        ft = ts[first[0]] if len(first) else None
        print(f"  run {r}: first change at t={ft}   count history: {comp}"
              f"{' ...' if len(runs) > 14 else ''}")
        rev = sum(1 for k in range(1, len(runs) - 1)
                  if runs[k][0] != 3 and runs[k + 1][0] == 3)
        print(f"          {len(runs)} plateaus, {rev} of them returned to 3 "
              f"(i.e. were blinks), final count {c[-1]}")
