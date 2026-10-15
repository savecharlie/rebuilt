#!/usr/bin/env python3
"""band.py -- the band is made of the one block form that belongs to every beat.
Does it stop growing, and at what width?

sweep32b.py found that under some implementations a charge-1 wall does not stay a
line: it widens into a region of blocks that carry no exclusive beat.

*** THE PREMISE OF THIS FILE WAS WRONG AND invade.py IS THE CORRECTION. ***
I wrote here that the region is "the form `100`, legal in all three alphabets".
It is not.  It is 100% the form `001`, which is legal in NO alphabet -- an
ILLEGAL phase that only a defect can manufacture.  `labels()` returns -1 for both
and I read every one as the all-beats kind.  The measurement below stands (the
region grows linearly and eats the ring, far above a measured chance-run null);
its interpretation is in invade.py.  Kept rather than deleted, because the null
control in it is the reason the growth is believable at all.

MEASURED, NOT SMOOTHED.  A band is defined on the RAW labels as a maximal run of
blocks carrying no exclusive beat.  The bulk is ~41% unlabelled, so runs occur by
chance: the null is measured on a wall-free locked ring of the same length under
the same triple, and the band must beat the null's longest run.  Without that
control any threshold would be a number I chose.
"""
import numpy as np

from tri import Kit

N, T = 600, 2500
BANDERS = [(36, 37, 13), (36, 37, 29), (38, 37, 13), (38, 37, 29)]
SHARP = [(36, 37, 5), (54, 37, 5)]


def longest_runs(lab):
    """-> array of (start, length) for every maximal run of unlabelled blocks,
    on the circle."""
    u = (lab < 0).astype(np.int8)
    if u.all():
        return np.array([[0, len(u)]])
    k = int(np.argmin(u))                       # a labelled block to cut at
    v = np.roll(u, -k)
    d = np.diff(np.r_[0, v, 0])
    s = np.flatnonzero(d == 1)
    e = np.flatnonzero(d == -1)
    return np.c_[(s + k) % len(u), e - s]


def null_longest(kit, seeds=8):
    """longest chance run of unlabelled blocks on a WALL-FREE ring, same size."""
    best = []
    for sd in range(seeds):
        s = kit.locked_ring(N, np.random.default_rng(770 + sd))
        m = 0
        for t in range(300):
            r = longest_runs(kit.labels(s, t))
            m = max(m, int(r[:, 1].max()) if len(r) else 0)
            s = kit.step(s)
        best.append(m)
    return max(best), float(np.mean(best))


if __name__ == "__main__":
    for tri in BANDERS + SHARP:
        kit = Kit(tri)
        nmax, nmean = null_longest(kit)
        widths = []
        for sd in range(4):
            s = kit.two_domain(N, 1, np.random.default_rng(31 + sd))
            w = []
            for t in range(T):
                r = longest_runs(kit.labels(s, t))
                w.append(int(r[:, 1].max()) if len(r) else 0)
                s = kit.step(s)
            widths.append(w)
        W = np.array(widths, float)
        late = W[:, -500:]
        tag = "BAND" if W[:, -1].mean() > nmax else "no band"
        print(f"\n=== {tri} ===")
        print(f"   null (no wall, same triple): longest chance run of unlabelled "
              f"blocks = {nmax} (mean over seeds {nmean:.1f})")
        qs = [int(W[:, min(t, T - 1)].mean()) for t in (50, 100, 250, 500, 1000, 2000, T - 1)]
        print(f"   widest run at t = 50/100/250/500/1000/2000/{T-1}: "
              + " / ".join(map(str, qs)))
        print(f"   last 500 steps: mean {late.mean():.1f}  sd across seeds "
              f"{late.mean(1).std(ddof=1):.1f}  per-seed finals "
              f"{[int(x) for x in W[:,-1]]}")
        print(f"   -> {tag}" + (f", saturated width {late.mean():.0f} blocks "
                                f"({late.mean()/nmax:.1f}x the chance run)" if tag == "BAND" else ""))
