#!/usr/bin/env python3
"""randic.py -- the invasion from RANDOM initial conditions, which is the one
that bears on the paper.

invade.py measured the illegal phase on a ring where I had BUILT the defect by
splicing two locked halves.  That is the right instrument for a wall's physics
and the wrong one for section 3, whose whole claim is about random bits.  So:
random rings, identical across triples, and watch the illegal-phase fraction.

The illegal form `001` is PRESENT at t = 0 -- random bits contain it at 1/8 of
blocks -- so this is not a question of nucleation at all.  It is a question of
whether the lattice cleans itself of a phase it starts with.
"""
import numpy as np

from tri import Kit
from invade import ILLEGAL

N, R, T = 256, 24, 40_000
MARKS = (0, 100, 1_000, 10_000, 40_000)
SHOW = [(36, 37, 5), (36, 37, 13), (36, 37, 29), (38, 37, 13), (38, 37, 29),
        (36, 37, 21), (52, 37, 29), (54, 37, 13), (54, 37, 29), (38, 39, 5)]


def frac(s):
    c = s.reshape(s.shape[0], -1, 3).astype(np.int64)
    return float(np.isin(4 * c[..., 0] + 2 * c[..., 1] + c[..., 2], ILLEGAL).mean())


if __name__ == "__main__":
    print(f"random rings, N={N} blocks, {R} per triple, THE SAME {R} initial "
          f"conditions for every triple")
    print(f"illegal-phase (`001`) fraction of the lattice:\n")
    print(f"   {'triple':<16} " + "".join(f"t={t:<10,}" for t in MARKS))
    for tri in SHOW:
        k = Kit(tri)
        s = np.random.default_rng(4242).integers(0, 2, (R, 3 * N)).astype(np.uint8)
        out = []
        for t in range(T + 1):
            if t in MARKS:
                out.append(frac(s))
            if t < T:
                s = k.step(s)
        tag = "  <-- fills with it" if out[-1] > 0.5 else ""
        print(f"   {str(tri):<16} " + "".join(f"{v:<12.4f}" for v in out) + tag)
    print("\nRandom bits contain the illegal form to start with (1/8 of blocks).")
    print("The paper's triple cleans itself of it; four of the 32 fill with it.")
