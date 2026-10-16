#!/usr/bin/env python3
"""chirality.py -- the (2,2,2)-never / (1,1,1)-always claim, measured over seeds.

NOTES.md's "best thing in this folder" says: three charge-(2,2,2) walls, 26 of 26
rings never resolved in 150,000 steps; three charge-(1,1,1) walls, 51 of 51
resolved, median 54,900. builttrap.py draws every run from ONE continuing
np.random.default_rng(31415) stream, so those are 26 and 51 particular lattices
and not a rate with an interval on it.

With independent per-seed lattices I found (2,2,2) rings locking at t = 500
(N=300) and t = 1400 (N=900), which "never in 150,000" does not allow. So the
rate gets measured.

The lock criterion is builttrap's own, which touches no labels: every block has
its middle and last bit zero (code 0 or 4). It is beat-sensitive -- a uniformly
locked ring is in that form at ONE of three beats -- so it is checked EVERY step,
which is what builttrap and race.py both do and what crosscheck.py's first
version failed to do.
"""
import json
import sys

import numpy as np

from wallcount import locked_ring, step, walls

W = 15


def ring_of_beats(N, beats, seed):
    rng = np.random.default_rng(seed)
    s = locked_ring(N, rng)
    k = len(beats)
    for i, off in enumerate(beats):
        b = locked_ring(N, rng)
        for _ in range(off):
            b = step(b)
        lo, hi = i * N // k, (i + 1) * N // k
        s[3 * lo:3 * hi] = b[3 * lo:3 * hi]
    return s


def lock_time(N, beats, seed, tmax):
    s = ring_of_beats(N, beats, seed)
    ch, n0 = walls(s, 0, W)
    q0 = sorted(int(c) for c in ch[ch != 0])
    for t in range(tmax):
        b = s.reshape(N, 3)
        if ((b[:, 1] == 0) & (b[:, 2] == 0)).all():
            return t, q0
        s = step(s)
    return None, q0


if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    tmax = int(sys.argv[2]) if len(sys.argv) > 2 else 150_000
    NS = int(sys.argv[3]) if len(sys.argv) > 3 else 40
    print(f"chirality.py -- N={N}, tmax={tmax:,}, {NS} independent seeds per cell")
    print("lock criterion is builttrap's own (no labels), checked EVERY step\n")
    out = {}
    for label, beats in (("(1,1,1)", (0, 1, 2)), ("(2,2,2)", (0, 2, 1))):
        times, bad = [], 0
        for sd in range(NS):
            t, q0 = lock_time(N, beats, 7000 + sd, tmax)
            if q0 != [int(label[1])] * 3:
                bad += 1
                continue
            times.append(t)
        res = [t for t in times if t is not None]
        print(f"  charge {label}: {len(res)}/{len(times)} rings locked within "
              f"{tmax:,} steps   (refused at t=0: {bad})")
        if res:
            print(f"            lock times: median {int(np.median(res)):,}  "
                  f"min {min(res):,}  max {max(res):,}")
        out[label] = dict(n=len(times), locked=len(res),
                          median=float(np.median(res)) if res else None,
                          times=[int(t) if t is not None else None for t in times])
    with open(f"chirality_N{N}.json", "w") as f:
        json.dump(out, f, indent=1)
    print(f"\nwrote chirality_N{N}.json")
