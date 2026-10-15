#!/usr/bin/env python3
"""pooled.py -- two triples are the SAME defect dynamics only if they agree on
EVERY wall configuration, not on one.

distinct.py found 24 classes on one charge-1 ring and 16 on one charge-2 ring.
A single configuration under-separates: it only exercises the free entries it
happens to reach.  Pool the partitions over many independent spliced rings (the
common refinement) and watch the class count climb.  The question is where it
stops: 32 would mean every one of the 32 implementations of the computation is
its own dynamical system once a defect exists.
"""
import hashlib

import numpy as np

from tri import MAIN_CUBE, Kit

N, T = 600, 150


def sig(kit, s0, T=T):
    h = hashlib.blake2b(digest_size=16)
    s = s0.copy()
    for _ in range(T):
        s = kit.step(s)
        h.update(s.tobytes())
    return h.digest()


if __name__ == "__main__":
    kits = {t: Kit(t) for t in MAIN_CUBE}
    base = kits[(36, 37, 5)]
    acc = {t: b"" for t in MAIN_CUBE}
    print(f"pooled over independent spliced rings (N={N} blocks, T={T} steps each)")
    print("  rings   charge-1 only   charge-2 only   both pooled")
    acc1 = {t: b"" for t in MAIN_CUBE}
    acc2 = {t: b"" for t in MAIN_CUBE}
    for i in range(1, 41):
        for charge, acc_c in ((1, acc1), (2, acc2)):
            s0 = base.two_domain(N, charge, np.random.default_rng(7000 + 97 * i + charge))
            for t, k in kits.items():
                acc_c[t] += sig(k, s0)
        for t in MAIN_CUBE:
            acc[t] = acc1[t] + acc2[t]
        n1 = len(set(acc1.values()))
        n2 = len(set(acc2.values()))
        nb = len(set(acc.values()))
        if i <= 10 or i % 5 == 0:
            print(f"  {i:>5}   {n1:>13}   {n2:>13}   {nb:>11}")
        if nb == 32:
            print(f"\n  -> all 32 separated after {i} ring(s) of each charge")
            break
    else:
        g = {}
        for t, v in acc.items():
            g.setdefault(v, []).append(t)
        print(f"\n  -> stalled at {len(set(acc.values()))} classes after 40 rings; "
              f"classes that never separated:")
        for c in g.values():
            if len(c) > 1:
                print(f"     {sorted(c)}")
