#!/usr/bin/env python3
"""front.py -- the invading illegal region has two edges. Is it still a WALL?

Every wall in this directory is classified by a Z3 charge: the jump in beat
phase across it, which must sum to 0 mod 3 around a ring.  The invading region
of `001` (invade.py) has no beat at all on the inside, so the usual question
"what is this wall's charge" has no obvious answer -- and yet the region has two
edges, and outside each edge is an ordinary locked domain with an ordinary beat.

THE QUESTION: is the phase difference across the WHOLE illegal region conserved
as the region grows?  If it is, the region is a fat wall -- the topological
charge survives while the structure that used to carry it does not.  If it is
not, the illegal phase is destroying the conservation law, which would be a
different and larger claim.

METHOD, and the trap it is built to avoid.  The phase labeller and the wall
detector are both built on the three alphabets and neither is meaningful inside
an illegal region, so NEITHER is used there.  The region is found on the RAW
codes (`001` is a code, not a label), and the phase is read only in the legal
lattice OUTSIDE it, at a fixed standoff, by majority over a window -- which is
the one place the labeller is valid.
"""
import numpy as np

from tri import Kit
from invade import ILLEGAL, codes

N, T, STANDOFF, WIN = 1800, 2200, 25, 40


def regions(ill):
    """maximal circular runs of True -> list of (start, length)."""
    if ill.all():
        return [(0, len(ill))]
    k = int(np.argmin(ill))
    v = np.roll(ill, -k).astype(np.int8)
    d = np.diff(np.r_[0, v, 0])
    s = np.flatnonzero(d == 1); e = np.flatnonzero(d == -1)
    return [((a + k) % len(ill), b - a) for a, b in zip(s, e)]


def phase_at(kit, s, t, j):
    """majority beat label over WIN blocks starting at j, legal lattice only."""
    lab = kit.labels(s, t)
    w = lab[(np.arange(WIN) + j) % N]
    w = w[w >= 0]
    if len(w) < WIN // 4:
        return -1
    return int(np.bincount(w).argmax())


if __name__ == "__main__":
    for tri in [(36, 37, 13), (38, 37, 29)]:
        kit = Kit(tri)
        print(f"\n=== {tri} ===")
        for sd in range(3):
            s = kit.two_domain(N, 1, np.random.default_rng(31 + sd))
            hist = []
            for t in range(T):
                if t and t % 200 == 0:
                    ill = np.isin(codes(s), ILLEGAL)
                    R = [r for r in regions(ill) if r[1] >= 20]
                    row = []
                    for (a, L) in R:
                        pl = phase_at(kit, s, t, (a - STANDOFF - WIN) % N)
                        pr = phase_at(kit, s, t, (a + L + STANDOFF) % N)
                        row.append((L, pl, pr, (pr - pl) % 3 if pl >= 0 and pr >= 0 else None))
                    hist.append((t, row))
                s = kit.step(s)
            print(f"  seed {sd}:")
            for t, row in hist[::2]:
                if not row:
                    print(f"    t={t:>5}  no illegal region >= 20 blocks")
                    continue
                print(f"    t={t:>5}  " + "   ".join(
                    f"[width {L:>4}  beat {pl}->{pr}  charge {c}]" for (L, pl, pr, c) in row))
