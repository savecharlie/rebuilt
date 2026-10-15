#!/usr/bin/env python3
"""bits.py -- the 32 implementations are the 32 assignments of 5 free bits.
Which bits does a wall actually feel, and how does the velocity depend on them?

The 5 free (residue, neighbourhood) entries are exactly the 5 bits that the
emulation does not force (reachable.py).  Writing each triple as that 5-bit word
turns "32 arbitrary triples" into a hypercube, and a velocity measured at each
vertex is a Boolean function on it.  Two things follow that a table of 32 rows
cannot show:

  - which single-bit flips move the velocity at all (the live bits), and
  - whether a bit is live only in the presence of another (masking), which
    would mean the wall's own structure decides which entries it ever reaches.
"""
import numpy as np

from tri import MAIN_CUBE, Kit
from freevisit import FREE

NAMES = [f"r{r}:{''.join(map(str,k))}" for r, k in FREE]
CODE = [(r, 4 * k[0] + 2 * k[1] + k[2]) for r, k in FREE]


def word(tri):
    """-> 5-bit tuple: the value of each free entry in this triple."""
    return tuple((tri[r] >> c) & 1 for r, c in CODE)


if __name__ == "__main__":
    print("the five free bits, and the rule each selects")
    for n, (r, c) in zip(NAMES, CODE):
        fam = {0: (36, 38, 52, 54), 1: (37, 39), 2: (5, 13, 21, 29)}[r]
        print(f"   {n:<12} residue {r}: " +
              "  ".join(f"{x}->{(x>>c)&1}" for x in fam))
    print()
    for t in MAIN_CUBE[:4] + MAIN_CUBE[-2:]:
        print(f"   {str(t):<16} -> {''.join(map(str, word(t)))}")
    assert len({word(t) for t in MAIN_CUBE}) == 32, "the 32 triples are not the 32 words"
    print("\n   the 32 triples are exactly the 32 five-bit words: ok")

    # ---- which free entries does each triple's own wall ever reach? --------
    print("\nfree entries VISITED at a charge-1 wall, per implementation")
    print("(a bit the wall never reaches cannot matter; this is measured, not reasoned)")
    N, T = 900, 400
    FM = np.zeros((3, 8), bool)
    for r, c in CODE:
        FM[r, c] = True
    res_all = np.tile(np.array([0, 1, 2]), N)
    out = {}
    for t_ in MAIN_CUBE:
        k = Kit(t_)
        seen = np.zeros(5, np.int64)
        for sd in range(3):
            s = k.two_domain(N, 1, np.random.default_rng(31 + sd))
            for _ in range(T):
                kk = (np.roll(s, 1, -1).astype(np.int64) << 2) | (s.astype(np.int64) << 1) \
                     | np.roll(s, -1, -1).astype(np.int64)
                for i, (r, c) in enumerate(CODE):
                    seen[i] += int(((res_all == r) & (kk == c)).sum())
                s = k.step(s)
        out[t_] = seen
    print(f"   {'triple':<16} {'word':<7} " + " ".join(f"{n:>12}" for n in NAMES))
    for t_ in MAIN_CUBE:
        w = "".join(map(str, word(t_)))
        print(f"   {str(t_):<16} {w:<7} " +
              " ".join(f"{v:>12,}" for v in out[t_]))
    np.save("freevisits_by_triple.npy",
            np.array([out[t] for t in MAIN_CUBE]))
