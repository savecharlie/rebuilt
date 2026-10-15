#!/usr/bin/env python3
"""gate256.py -- run the GOOD gate at N = 256, on the run that was gated by the bad one.

lock32.json was produced by a process that was already in memory when I replaced
the median-reproduction gate (law 30) with the instrument-against-instrument one.
It passed the gate it had, which does not mean much. The MEASUREMENT code is the
same; only the gate differs. So rather than redo a two-hour sweep, redo the GATE:
my locking criterion (every block in one of the three alphabets) against
lockon.py's narrower `[x,0,0]`, on the SAME 160 rings, at N = 256.

Tolerance derived, as in lock32_N128.py: 11 steps of measured early-firing (my
criterion is not absorbing) plus CHECK of quantisation.
"""
import json

import numpy as np

from tri import Kit, MAIN_CUBE
from lock32 import N, RUNS, TMAX, CHECK, locked

kb = Kit(MAIN_CUBE[0])
rng = np.random.default_rng(5150)
ic = rng.integers(0, 2, (RUNS, 3 * N)).astype(np.uint8)


def strict(x):
    b = x.reshape(x.shape[0], -1, 3)
    return ((b[:, :, 1] == 0) & (b[:, :, 2] == 0)).all(1)


def sweep(crit, three_step):
    s = ic.copy(); row = np.arange(RUNS); tl = np.full(RUNS, -1); t = 0
    while len(row) and t < TMAX:
        if three_step:
            for _ in range(3):
                s = kb.step(s); t += 1
                g = crit(s)
                if g.any():
                    tl[row[g]] = t; s = s[~g]; row = row[~g]
                if not len(row):
                    break
            if not len(row):
                break
            for _ in range(CHECK - 3):
                s = kb.step(s); t += 1
        else:
            for _ in range(CHECK):
                s = kb.step(s); t += 1
            g = crit(s)
            if g.any():
                tl[row[g]] = t; s = s[~g]; row = row[~g]
    return tl


if __name__ == "__main__":
    print(f"GATE L at N={N}, {RUNS} rings, budget {TMAX:,}")
    mine = sweep(locked, False)
    theirs = sweep(strict, True)
    agree = ((mine >= 0) == (theirs >= 0)).all()
    both = (mine >= 0) & (theirs >= 0)
    lag = np.abs(theirs[both] - mine[both])
    TOL = 11 + CHECK
    ok = agree and lag.max() <= TOL
    print(f"   same locked set: {agree}   ({(mine>=0).sum()} vs {(theirs>=0).sum()} "
          f"of {RUNS})")
    print(f"   max |t_lock difference|: {lag.max()}   tolerance {TOL} "
          f"(11 early-firing + {CHECK} quantisation)")
    print(f"   -> {'PASS' if ok else 'FAIL'}")
    print(f"   fraction locked: mine {(mine>=0).mean():.3f}, strict "
          f"{(theirs>=0).mean():.3f}   (lock32.json reports "
          f"{json.load(open('lock32.json'))['data']['(36, 37, 5)']['frac']:.3f})")
