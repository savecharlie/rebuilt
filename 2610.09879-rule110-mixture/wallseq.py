#!/usr/bin/env python3
"""wallseq.py -- the walls and their charges, read straight off the raw labels.

THE CORRECTION THIS FILE IS. burnin.py asked "is the shortest plateau at the
intermediate phase", which was built for a short strip hiding inside a wall.
Once a third domain GROWS, the shortest plateau is some other domain and the
criterion silently stops firing -- it reported 19% while logging 2,538
"off-phase" frames that were the same event seen through a filter whose
membership depends on the answer. Eighth instrument on this paper to lie, and
a repeat species: a SELECTING filter (fire 343's first five were the same
shape -- sorted, minimum, surviving, readable, slot-zero).

What is measured instead, and it needs no criterion at all: the cyclic sequence
of unambiguous phase plateaus around the ring gives the walls and their charges
EXACTLY. Plateau phases (p0, p1, ... pk) -> charges ((p1-p0) %3, ..., (p0-pk) %3).
Two plateaus is two walls; three is three. The charges must sum to 0 mod 3, and
that sum is a free CHECK on every frame -- if it ever fails, the reading is
wrong and the frame is refused rather than counted.

The physics question: a ring built with two walls (+1, +2) is observed by
fire 340 to end with THREE walls all of charge +1 -- the paper's trapped state.
2 = 1 + 1, so that transition IS the charge-2 wall splitting. Does it happen,
how often, and does the ring ever go back?
"""
import json

import numpy as np

from split import plateaus
from tri import Kit

BASE = (36, 37, 5)
N, BURN = 600, 1_000


def wall_charges(kit, s, t):
    """(charges, phases) read off the plateau sequence. None if the sum fails."""
    pl = plateaus(kit, s, t)
    if len(pl) < 2:
        return (), tuple(p for p, _, _ in pl)
    ph = [p for p, _, _ in pl]
    ch = tuple((ph[(i + 1) % len(ph)] - ph[i]) % 3 for i in range(len(ph)))
    if sum(ch) % 3 != 0:
        return None, tuple(ph)          # refused: the reading is inconsistent
    return ch, tuple(ph)


def arc(kit, d, ahead, seed):
    rng = np.random.default_rng(seed)
    s = kit.locked_ring(N, rng)
    b = kit.locked_ring(N, np.random.default_rng(seed))
    for _ in range(ahead):
        b = kit.step(b)
    s[3 * 200:3 * (200 + d)] = b[3 * 200:3 * (200 + d)]
    return s


def trace(kit, s, T=30_000, every=10):
    hist, refused, first3, last2 = {}, 0, None, None
    for _ in range(BURN):
        s = kit.step(s)
    for t in range(BURN, BURN + T):
        if t % every == 0:
            ch, _ = wall_charges(kit, s, t)
            if ch is None:
                refused += 1
            else:
                key = tuple(sorted(ch))
                hist[key] = hist.get(key, 0) + 1
                if len(ch) >= 3 and first3 is None:
                    first3 = t
                if len(ch) == 2:
                    last2 = t
        s = kit.step(s)
    return hist, refused, first3, last2


if __name__ == "__main__":
    kit = Kit(BASE)
    print("wallseq.py -- wall charges read exactly off the plateau sequence")
    print(f"triple {BASE}, N={N}, burn-in {BURN:,}, 30,000 steps every 10")
    print("charges must sum to 0 mod 3; a frame that fails is REFUSED, not counted\n")
    out = {}
    for ahead in (1, 2):
        for d in (200, 300):
            print(f"=== arc d={d} blocks, {ahead} beat(s) ahead "
                  f"(built walls: +{ahead}, +{(3-ahead)%3}) ===")
            for sd in range(8):
                hist, refused, f3, l2 = trace(kit, arc(kit, d, ahead, 300 + sd))
                tot = sum(hist.values())
                desc = "  ".join(
                    f"{'+'.join(map(str, k))}:{100*v/tot:.0f}%"
                    for k, v in sorted(hist.items(), key=lambda kv: -kv[1])[:4])
                n3 = sum(v for k, v in hist.items() if len(k) >= 3)
                print(f"   seed {sd}: {100*n3/tot:5.1f}% of frames have >=3 walls"
                      f"   refused {refused:4d}   first3 "
                      f"{f3 if f3 else '-':>6}   last2 {l2 if l2 else '-':>6}")
                print(f"            {desc}")
                out[f"d{d}_a{ahead}_s{sd}"] = dict(
                    frac3=n3 / tot, refused=refused, first3=f3, last2=l2,
                    hist={'+'.join(map(str, k)): v for k, v in hist.items()})
            print()
    with open("wallseq.json", "w") as f:
        json.dump(out, f, indent=1)
    print("wrote wallseq.json")
