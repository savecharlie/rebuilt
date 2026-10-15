#!/usr/bin/env python3
"""p17.py -- the registered question from fire 342: does the INTEGRATED core
separate by charge, with more seeds than the five that left it ambiguous?

core.py, five seeds: peak 0.085 +- 0.007 (charge 1) vs 0.058 +- 0.003 (charge 2),
3.5 sigma; integral 0.193 +- 0.036 vs 0.161 +- 0.006, not separated, because one
charge-1 seed of five gave nearly double the others.

Twelve independent seeds, and the per-seed numbers are printed, because a heavy
tail is a thing you look at rather than summarise.  Same identity-tracked
instrument as sweep32b.py (modal charge and modal domain phases over the run),
which is NOT core.py's per-frame instrument -- and the two disagree on the
absolute peak by 30%, so that disagreement is reported here rather than hidden.
"""
import numpy as np

from tri import Kit
from sweep32b import one

if __name__ == "__main__":
    kit = Kit((36, 37, 5))
    res = {1: [], 2: []}
    for sd in range(12):
        r = one(kit, sd)
        if r is None:
            print(f"   seed {sd}: dropped")
            continue
        for c, (v, mid, amb) in r.items():
            res[c].append((sd, float(mid.max()), float(mid.sum()), float(v)))
    for c in (1, 2):
        a = res[c]
        print(f"\n=== charge {c}: {len(a)} seeds ===")
        for sd, pk, ig, v in a:
            print(f"   seed {sd:>2}  peak {pk:.4f}   integral {ig:.4f}   v {v:+.4f}")
        pk = np.array([x[1] for x in a]); ig = np.array([x[2] for x in a])
        print(f"   peak     {pk.mean():.4f} +- {pk.std(ddof=1)/np.sqrt(len(pk)):.4f}"
              f"   (sd {pk.std(ddof=1):.4f}, max/median {pk.max()/np.median(pk):.2f})")
        print(f"   integral {ig.mean():.4f} +- {ig.std(ddof=1)/np.sqrt(len(ig)):.4f}"
              f"   (sd {ig.std(ddof=1):.4f}, max/median {ig.max()/np.median(ig):.2f})")
        res[c] = (pk, ig)
    for name, i in (("peak", 0), ("integral", 1)):
        a, b = res[1][i], res[2][i]
        z = (a.mean() - b.mean()) / np.hypot(a.std(ddof=1) / np.sqrt(len(a)),
                                             b.std(ddof=1) / np.sqrt(len(b)))
        print(f"\n{name:>8}: charge1 {a.mean():.4f} vs charge2 {b.mean():.4f}  "
              f"-> z = {z:+.2f}  {'SEPARATED' if abs(z) > 2 else 'not separated'}")
