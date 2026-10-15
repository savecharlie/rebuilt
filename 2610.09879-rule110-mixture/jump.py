#!/usr/bin/env python3
"""jump.py -- how far a DETECTION moves between samples, measured not assumed.

tracked.Track capped the per-sample jump at MAXV*every*3+2 blocks, built from the
wall's measured velocity ceiling of 0.25 blocks/step. That is the wall's motion.
It is not the DETECTOR's: the position comes from a majority vote over a 15-block
box on a label field that is 40% holes, so it scatters by blocks even when the
wall is standing still. With every=5 the cap was 5.75 and real jumps of 6 were
throwing the tracker off a perfectly healthy ring -- at t=2100 the detections ran
[86,144,239] -> [86,144,233] and the match refused it.

Walls here sit 60+ blocks apart, so nearest-neighbour matching between two
3-detection frames is unambiguous no matter how big the jump. Measure it.
"""
import numpy as np
from wallcount import step, walls
from tracked import _ring3

N, W = 300, 15
for every in (1, 5, 25):
    allj = []
    for offs, name in (((0, 1, 2), "(1,1,1)"), ((0, 2, 1), "(2,2,2)")):
        rng = np.random.default_rng(342342)
        js = []
        for r in range(4):
            s = _ring3(N, [0, 60, 150], offs, rng)
            prev = None
            for t in range(4001):
                if t % every == 0:
                    ch, _ = walls(s, t, W)
                    p = np.flatnonzero(ch)
                    if len(p) == 3:
                        if prev is not None:
                            d = np.abs(p[:, None] - prev[None, :])
                            d = np.minimum(d, N - d)
                            js.extend(d.min(1).tolist())
                        prev = p
                s = step(s)
        js = np.array(js)
        allj.extend(js.tolist())
        print(f"every={every:>2} {name}: n={len(js):6d}  mean {js.mean():5.2f}  "
              f"p50 {np.percentile(js,50):4.0f}  p99 {np.percentile(js,99):4.0f}  "
              f"p99.9 {np.percentile(js,99.9):5.0f}  max {js.max():4d}")
    a = np.array(allj)
    print(f"          pooled p99.9 = {np.percentile(a,99.9):.0f}, max = {a.max()}, "
          f"wall motion over {every} steps = {0.25*every:.2f}\n")
