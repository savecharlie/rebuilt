#!/usr/bin/env python3
"""nullwalk.py -- what SORTED gap means do when there is provably no force.

Three walls on a ring of N = 300 blocks. Each does an independent random walk.
There is no interaction of any kind: I wrote the generator, the three increments
are drawn independently, and nothing in the code can see another wall. So the
ground truth is known -- zero force -- and any "relaxation toward equal spacing"
that comes out of this is an artifact of how the gaps were summarised.

This is the pre-test for P11. Fire 341 reported a charge-(2,2,2) triple started
at gaps (60, 90, 150) relaxing to (79, 97, 124) inside 8,000 steps and holding,
and called it a restoring force. The number was a mean of np.sort()ed gaps over
the runs that survived. The question here is what that same summary does to a
forceless walk.

Absorbing: a run ends when any gap reaches 0 (two walls have met). Reported both
ways -- over survivors only, like repel.py did, and over all runs with the
survivor count beside it.
"""
import numpy as np

N = 300
G0 = (60.0, 90.0, 150.0)
RUNS = 2000
T = 40000
EVERY = 4000


def walk(sigma, runs=RUNS, T=T, every=EVERY, seed=0):
    """-> t, gaps[run, frame, 3] by IDENTITY (NaN after absorption), alive[frame]."""
    rng = np.random.default_rng(seed)
    # wall positions: x0 = 0, x1 = 60, x2 = 150 on the ring; gaps are differences
    x = np.zeros((runs, 3))
    x[:, 1], x[:, 2] = G0[0], G0[0] + G0[1]
    nf = T // every + 1
    out = np.full((runs, nf, 3), np.nan)
    dead = np.zeros(runs, bool)

    def gaps(x):
        return np.stack([x[:, 1] - x[:, 0], x[:, 2] - x[:, 1],
                         N + x[:, 0] - x[:, 2]], -1)

    out[:, 0, :] = gaps(x)
    for t in range(1, T + 1):
        x += rng.normal(0.0, sigma, x.shape)
        g = gaps(x)
        dead |= (g <= 0).any(-1)
        if t % every == 0:
            k = t // every
            out[:, k, :] = np.where(dead[:, None], np.nan, g)
    return np.arange(0, T + 1, every), out, (~np.isnan(out[:, :, 0])).sum(0)


if __name__ == "__main__":
    print(f"three independent walks, N={N}, gaps built at {G0}, no force by construction\n")
    for sigma in (0.05, 0.1, 0.2, 0.4):
        t, g, alive = walk(sigma)
        srt = np.sort(g, -1)
        print(f"sigma = {sigma}  (rms single-gap displacement over 8000 steps "
              f"= {sigma*np.sqrt(2*8000):.0f} blocks)")
        print(f"  {'t':>7} {'alive':>6}   sorted means (survivors)        "
              f"identity means (survivors)      spread srt  id")
        for k in (0, 1, 2, 5, 10):
            if k >= len(t):
                continue
            m = np.nanmean(srt[:, k, :], 0)
            mi = np.nanmean(g[:, k, :], 0)
            print(f"  {t[k]:>7} {alive[k]:>6}   "
                  f"({m[0]:5.1f},{m[1]:5.1f},{m[2]:5.1f})   "
                  f"({mi[0]:5.1f},{mi[1]:5.1f},{mi[2]:5.1f})   "
                  f"{m[2]-m[0]:9.1f} {mi.max()-mi.min():4.1f}")
        print()
