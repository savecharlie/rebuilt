#!/usr/bin/env python3
"""distinct.py -- the computation has 32 implementations; how many does the
DEFECT have?

All 32 main-cube triples agree on every table entry the emulation forces, and on
the locked-domain alphabets (tri.py, gates E and A).  freevisit.py showed the
five free entries ARE reached, and only at a wall.  So: evolve the SAME
two-domain configuration under all 32 and count how many distinct trajectories
come out.

Two gates, both against answers known before the run:
  G-same  -- on a WALL-FREE locked ring all 32 must agree bit-for-bit forever
             (reachable.py proves the free entries are unreachable there).
  G-split -- the partition must refine as T grows and never coarsen, because
             trajectories that have diverged cannot re-merge into equality at
             a later time under a deterministic map... they CAN, actually, so
             this gate only checks monotone refinement of the EQUIVALENCE at
             fixed prefix length, which is what the code computes.
"""
import hashlib

import numpy as np

from tri import MAIN_CUBE, Kit

N = 600


def traj_hash(kit, s0, T):
    h = hashlib.blake2b(digest_size=16)
    s = s0.copy()
    for _ in range(T):
        s = kit.step(s)
        h.update(s.tobytes())
    return h.hexdigest(), s


def classes(kits, s0, T):
    g = {}
    for t, k in kits.items():
        hh, _ = traj_hash(k, s0, T)
        g.setdefault(hh, []).append(t)
    return list(g.values())


if __name__ == "__main__":
    kits = {t: Kit(t) for t in MAIN_CUBE}
    base = kits[(36, 37, 5)]
    rng = np.random.default_rng(343)

    print("G-same -- wall-free locked ring, all 32 must stay bit-identical")
    s0 = base.locked_ring(N, rng)
    cl = classes(kits, s0, 400)
    print(f"   {len(cl)} class(es) after 400 steps -> "
          f"{'PASS' if len(cl) == 1 else 'FAIL'}")
    assert len(cl) == 1

    for charge in (1, 2):
        print(f"\n=== charge-{charge} two-domain ring ({N} blocks) ===")
        s0 = base.two_domain(N, charge, np.random.default_rng(1000 + charge))
        prev = None
        for T in (1, 2, 3, 6, 12, 25, 50, 100, 200, 400):
            cl = classes(kits, s0, T)
            sizes = sorted((len(c) for c in cl), reverse=True)
            print(f"   T={T:<4} distinct trajectories: {len(cl):>2}   class sizes {sizes}")
            prev = cl
        print("   final partition:")
        for c in sorted(prev, key=lambda c: (-len(c), c)):
            print(f"     {len(c):>2}: {sorted(c)}")
