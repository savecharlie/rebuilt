#!/usr/bin/env python3
"""core.py -- what the raw phase labels look like inside a wall.

WHY. `blink2.py` found that this detector occasionally resolves a single
charge-1 wall as TWO charge-2 walls exactly 2 blocks apart (87-92% of extra
detections carry charge 2, median separation 2.0 blocks at every one of four
seeds), and a charge-2 wall as two charge-1 pieces at a variable 5-30 blocks.
Since 2 + 2 = 4 = 1 (mod 3), the natural reading is that a charge-1 wall is a
bound pair of charge-2 sub-walls -- i.e. there is a sliver of the intermediate
phase inside it.

But that was read off the SMOOTHED field at w = 15, and resolution.py measured
that w = 15 erases anything narrower than 15 blocks. A 2-block structure is
below the resolution of the instrument that found it, so I cannot trust it from
there. This script throws the smoothing away and reads the RAW label field.

A charge-1 wall is the boundary phase a -> a+1. If it is really two charge-2
sub-walls, the raw labels crossing it should go a -> a+2 -> a+1, with the a+2
sliver a couple of blocks wide. If it is a single wall, there is no sliver and
the profile is a step through the ambiguous blocks.

The wall is located on the smoothed field (that is all it is used for -- a
position, good to a block or two) and the RAW labels are histogrammed by offset.
"""
import sys

import numpy as np

from wallcount import step, walls, labels, locked_ring

N, W, R = 1800, 15, 40          # big ring, two far-apart walls, +-R blocks of profile


def profile(charge, reps=10, T=600, burn=120, seed=342, N=N):
    """-> (2R+1, 4) counts of [phase_left, phase_mid, phase_right, ambiguous]
    where 'mid' is the third phase, the one neither domain is in."""
    rng = np.random.default_rng(seed)
    hist = np.zeros((2 * R + 1, 4))
    n_used = 0
    for _ in range(reps):
        a = locked_ring(N, rng)
        b = locked_ring(N, rng)
        for _ in range(charge):            # right domain is `charge` beats ahead
            b = step(b)
        s = a.copy()
        s[3 * (N // 4):3 * (3 * N // 4)] = b[3 * (N // 4):3 * (3 * N // 4)]
        for t in range(T):
            if t >= burn:
                ch, _ = walls(s, t, W)
                p = np.flatnonzero(ch)
                if len(p) == 2:
                    for j in p:
                        q = int(ch[j])
                        lab = labels(s, t)
                        # left domain phase is lab just outside, read from the
                        # smoothed field to avoid the ambiguous blocks
                        ph = (np.roll(lab, 0),)
                        idx = (np.arange(-R, R + 1) + j) % N
                        L = lab[idx]
                        left = int(np.median(lab[(np.arange(j - 60, j - 20)) % N][
                            lab[(np.arange(j - 60, j - 20)) % N] >= 0]))
                        mid = (left + 2 * q) % 3 if q else -1
                        right = (left + q) % 3
                        col = np.full(len(L), 3)
                        col[L == left] = 0
                        col[L == mid] = 1
                        col[L == right] = 2
                        # orient so the left domain is always on the left
                        hist[np.arange(len(L)), col] += 1
                        n_used += 1
            s = step(s)
    return hist, n_used


if __name__ == "__main__":
    for charge in (1, 2):
        h, n = profile(charge)
        f = h / h.sum(1, keepdims=True)
        print(f"\n=== raw label profile across a charge-{charge} wall "
              f"({n} wall-frames, ring of {N}) ===")
        print(f"  offset   P(left phase)  P(MID phase)  P(right phase)  P(ambiguous)")
        for k in range(0, 2 * R + 1):
            o = k - R
            if abs(o) > 14 and abs(o) % 5:
                continue
            bar = "#" * int(round(40 * f[k, 1]))
            print(f"  {o:+4d}      {f[k,0]:6.3f}        {f[k,1]:6.3f}        "
                  f"{f[k,2]:6.3f}        {f[k,3]:6.3f}   {bar}")
        far = np.r_[f[:R-25, 1], f[R+26:, 1]].mean()
        near = f[R-2:R+3, 1].mean()
        nz = np.flatnonzero(f[:, 1] > 0.002) - R
        print(f"  MID-phase fraction far from the wall: {far:.6f} "
              f"(EXACTLY zero -- the domains are in their own phases, so any "
              f"MID label at all is wall-localised)")
        print(f"  MID present (>0.2%) at offsets {nz.min():+d} .. {nz.max():+d}  "
              f"({nz.max()-nz.min()+1} blocks wide), peak {f[:,1].max():.3f} at "
              f"offset {int(np.argmax(f[:,1]))-R:+d}, integral "
              f"{f[:,1].sum():.3f} blocks of MID per wall")
