#!/usr/bin/env python3
"""resolution.py -- what is the smallest domain the wall counter can see?

G1 in wallcount.py FAILED on narrow arcs and that is the gate working, not a
bug: smoothing over w blocks cannot resolve a domain narrower than ~w, and too
small a w invents walls out of the 41% unlabelled noise. So measure the limit
instead of tuning around it. For each w, the smallest arc width D at which the
count is EXACT on every one of `reps` random rings.
"""
import numpy as np
from wallcount import walls, locked_ring, step, encode


def arc_ring(N, width, off, rng, narcs=1):
    s = locked_ring(N, rng)
    slot = N // narcs
    for i in range(narcs):
        b = locked_ring(N, rng)
        for _ in range(off):
            b = step(b)
        lo = i * slot + (slot - width) // 2
        s[3 * lo:3 * (lo + width)] = b[3 * lo:3 * (lo + width)]
    return s


rng = np.random.default_rng(7)
N = 600
reps = 40
print(f"ring N={N} blocks, {reps} random rings per cell, arcs evenly spaced")
print("exact count required on ALL reps. '.' = exact, 'x' = wrong somewhere\n")
widths = [3, 4, 5, 6, 8, 10, 13, 16, 20, 26, 33, 42, 54, 70]
for narcs in (1, 3):
    want = 2 * narcs
    print(f"--- {narcs} arc(s), want {want} walls ---")
    print("  w \\ D " + "".join(f"{d:>5}" for d in widths))
    for w in (5, 7, 9, 11, 15, 21, 31):
        row = ""
        for D in widths:
            if narcs * D > N - narcs * 4:
                row += "    -"; continue
            good = True
            for r in range(reps):
                for off in (1, 2):
                    s = arc_ring(N, D, off, rng, narcs)
                    if int(walls(s, 0, w)[1]) != want:
                        good = False; break
                if not good:
                    break
            row += "    ." if good else "    x"
        print(f"  {w:>3}   {row}")
    print()
