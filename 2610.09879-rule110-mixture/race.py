#!/usr/bin/env python3
"""race.py -- is the invasion a RACE, and can annihilation win it?

front.py: under the four invading triples, a wall nucleates an illegal region that
grows at a fixed rate and conserves the charge it was born with, and two such
regions do not annihilate when they meet -- they merge, and the ring ends
uniformly illegal.

But a thin charge-1 and a thin charge-2 wall DO annihilate; that is the whole of
section 3's mechanism. So there ought to be a race: a pair that meets while still
thin should cancel and leave clean lattice, and a pair that widens first should
not. If that is right there is a critical SEPARATION below which an invaded ring
still locks, and lock32's zeros would be a statement about typical wall spacing
rather than about the triple.

If it is NOT right -- if even a pair built two blocks apart ends illegal -- then
the invasion is unconditional and a single defect anywhere is terminal.

CONSTRUCTION: an arc of `d` blocks one beat ahead inside an otherwise locked
ring, which is exactly two walls (charge 1 and charge 2) at separation d.
Outcome is read off the RAW codes, never the labeller: `locked` = no illegal
block and every block a codeword; `illegal` = any illegal block left at the end.
"""
import numpy as np

from tri import Kit
from invade import ILLEGAL, codes

N, T = 900, 20_000
SEPS = [2, 3, 5, 8, 12, 20, 32, 50, 80, 130, 200]


def arc(kit, N, d, rng, ahead=1):
    s = kit.locked_ring(N, rng)
    b = kit.locked_ring(N, rng)
    for _ in range(ahead):
        b = kit.step(b)
    lo = N // 3
    s[3 * lo:3 * (lo + d)] = b[3 * lo:3 * (lo + d)]
    return s


def outcome(kit, s, T=T):
    for t in range(T):
        c = codes(s)
        ill = np.isin(c, ILLEGAL)
        if not ill.any() and np.isin(c, [0, 4]).all():
            return "LOCKED", t
        s = kit.step(s)
    c = codes(s)
    f = float(np.isin(c, ILLEGAL).mean())
    return (f"illegal {f:.2f}" if f > 0.5 else f"open (illegal {f:.2f})"), T


if __name__ == "__main__":
    print(f"one arc of d blocks, one beat ahead = exactly two walls at separation d")
    print(f"ring {N} blocks, {T:,} steps, 6 seeds per cell\n")
    print(f"   {'triple':<14} " + "".join(f"d={d:<8}" for d in SEPS))
    for tri in [(36, 37, 5), (36, 37, 13), (38, 37, 29), (52, 37, 29)]:
        kit = Kit(tri)
        cells = []
        for d in SEPS:
            n_lock = 0
            for sd in range(6):
                o, _ = outcome(kit, arc(kit, N, d, np.random.default_rng(900 + sd)))
                n_lock += (o == "LOCKED")
            cells.append(f"{n_lock}/6")
        print(f"   {str(tri):<14} " + "".join(f"{c:<10}" for c in cells))
    print("\ncells are 'rings that reached a legal codeword configuration / rings tried'")
