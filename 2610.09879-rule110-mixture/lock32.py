#!/usr/bin/env python3
"""lock32.py -- P18: does the paper's self-organisation result survive a change
of implementation?

Section 3 of arXiv:2610.09879: random initial conditions LOCK ON -- reach a
configuration that is a legal codeword at some phase -- 91-100% of the time for
N >= 256.  That is what makes the construction a computer that assembles itself
rather than one that has to be prepared.  lockon.py reproduced it for (36,37,5):
0.933 +- 0.023 at N = 256, median t_lock 13,295.

Locking requires every wall to leave the ring.  sweep32b.py measured that under
many implementations a wall PINS: velocity 0.0008 +- 0.0005 against the
baseline's -0.245, either as a band of the ambiguous form `100` or as a
persistent sliver of the third beat.

ALL 32 IMPLEMENTATIONS ARE STEPPED IN ONE ARRAY, on the SAME 160 initial
conditions.  Run one at a time this takes four hours, because every ring of a
non-locking implementation burns the whole budget; batched it is one numpy call
per step for all of them.  Sharing the initial conditions also removes the
seed-to-seed term from every pairwise comparison -- each implementation is asked
the same 160 questions.

BUDGET.  Fixed at 200,000 steps, about 15x the measured median.  The gate below
reports that the baseline locks 0.863 at this budget against 0.933 at 16 N^2, so
the budget IS binding at the 7% level and every fraction here is "locked within
200,000 steps".  What does not depend on the budget is the median lock time of
the runs that do lock, and that is what GATE L tests.

"Locked" needs no tuning: every three-site block lies in ONE of the three
alphabets.  The state is absorbing (the emulation is exact for all 32, tri.py
gate E), so testing every CHECK steps coarsens t_lock and cannot change the
fraction.
"""
import json
import sys
import time

import numpy as np

from tri import MAIN_CUBE, Kit, BASE_ALPH

N, RUNS, TMAX, CHECK = 256, 160, 200_000, 25
ALPHSET = [np.array(sorted(4 * x + 2 * y + z for (x, y, z) in A))
           for A in BASE_ALPH.values()]
TABS = np.stack([Kit(t, check=False).TAB for t in MAIN_CUBE])      # (32, 3, 8)
RES = np.tile(np.array([0, 1, 2], np.int64), N)                    # (3N,)


def step_many(s, tid):
    k = (np.roll(s, 1, -1).astype(np.int64) << 2) | (s.astype(np.int64) << 1) \
        | np.roll(s, -1, -1).astype(np.int64)
    return TABS[tid[:, None], RES[None, :], k]


def locked(s):
    b = s.reshape(s.shape[0], -1, 3).astype(np.int64)
    code = 4 * b[..., 0] + 2 * b[..., 1] + b[..., 2]
    ok = np.zeros(s.shape[0], bool)
    for A in ALPHSET:
        ok |= np.isin(code, A).all(1)
    return ok


if __name__ == "__main__":
    nt = len(MAIN_CUBE)
    rng = np.random.default_rng(5150)
    ic = rng.integers(0, 2, (RUNS, 3 * N)).astype(np.uint8)        # SHARED
    s = np.repeat(ic[None], nt, 0).reshape(nt * RUNS, 3 * N)
    tid = np.repeat(np.arange(nt), RUNS)
    row = np.arange(nt * RUNS)
    tlock = np.full(nt * RUNS, -1)
    t = 0
    t0 = time.time()
    while len(row) and t < TMAX:
        for _ in range(CHECK):
            s = step_many(s, tid)
            t += 1
        m = locked(s)
        if m.any():
            tlock[row[m]] = t
            s = s[~m]; tid = tid[~m]; row = row[~m]
        if t % 5000 == 0:
            print(f"   t={t:>7,}  {len(row):>6,} of {nt*RUNS:,} rings still open"
                  f"   ({time.time()-t0:.0f}s)", flush=True)

    tl = tlock.reshape(nt, RUNS)
    f0 = (tl[0] >= 0).mean()
    med0 = int(np.median(tl[0][tl[0] >= 0]))
    g = abs(med0 - 13295) / 13295 < 0.20
    print(f"\nGATE L -- the baseline's median lock time, which does not depend "
          f"on the budget")
    print(f"   (36,37,5): median t_lock {med0:,}   ref 13,295   -> "
          f"{'PASS' if g else 'FAIL'}")
    print(f"   fraction locked within {TMAX:,} steps: {f0:.3f}   "
          f"(lockon.py got 0.933 +- 0.023 at a budget 5x longer: the budget IS "
          f"binding at the 7% level and every number below is at THIS budget)\n")
    if not g:
        sys.exit(1)
    out = {}
    for i, tri in enumerate(MAIN_CUBE):
        a = tl[i]
        f = (a >= 0).mean()
        med = int(np.median(a[a >= 0])) if (a >= 0).any() else -1
        se = float(np.sqrt(f * (1 - f) / RUNS))
        out[str(tri)] = {"frac": float(f), "se": se, "median_tlock": med,
                         "n_locked": int((a >= 0).sum())}
        print(f"   {str(tri):<16} locked {f:.3f} +- {se:.3f}   "
              f"median t_lock {med if med >= 0 else 'none':>8}")
    json.dump({"N": N, "runs": RUNS, "tmax": TMAX, "shared_ic": True,
               "data": out}, open("lock32.json", "w"))
    print("\nwrote lock32.json")
