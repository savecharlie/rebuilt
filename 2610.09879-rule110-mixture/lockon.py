#!/usr/bin/env python3
"""lockon.py -- section 3 of arXiv:2610.09879: do RANDOM initial conditions
self-organise into rule 110?

The paper: for N >= 256 blocks, ~91-100% of runs lock into rule 110, locking time
~ N^2, and 0-9% are trapped -- because the lattice breaks into regions emulating
rule 110 in different PHASES, and the walls between them move ballistically; a
left-mover pins a domain and takes over, a set of right-movers never meet.

"Locked" here = every block of three sites is a legal codeword, [1,0,0] or
[0,0,0]. Once that holds it holds forever at that phase, because the emulation is
exact (verify.py). That is the operational definition and it needs no tuning.

Runs are batched and retired as they lock, so the cost is set by the stragglers.
"""
import sys
import time

import numpy as np

RULES = (36, 37, 5)
TAB = np.array([[(r >> k) & 1 for k in range(8)] for r in RULES], np.uint8)


def step(s, res):
    L = s.shape[1]
    k = (np.roll(s, 1, 1).astype(np.uint8) << 2) | (s << 1) | np.roll(s, -1, 1)
    return TAB[res, k]


def run(N, nruns, tmax, seed=0):
    rng = np.random.default_rng(seed)
    L = 3 * N
    res = np.tile(np.array([0, 1, 2], np.int64), N)
    s = rng.integers(0, 2, (nruns, L), dtype=np.uint8)
    alive = np.arange(nruns)
    lock_t = np.full(nruns, -1, np.int64)
    t0 = time.time()
    for t in range(tmax + 1):
        b = s.reshape(len(s), N, 3)
        legal = (b[:, :, 1] == 0) & (b[:, :, 2] == 0)       # [x,0,0]
        done = legal.all(1)
        if done.any():
            lock_t[alive[done]] = t
            s = s[~done]
            alive = alive[~done]
            if len(s) == 0:
                break
        s = step(s, res)
    locked = lock_t >= 0
    return lock_t, locked, time.time() - t0


print(f"{'N':>6} {'runs':>5} {'tmax':>8} {'locked':>8} {'frac':>7} "
      f"{'median t_lock':>14} {'t_lock/N^2':>11} {'wall s':>7}")
for N, nruns, mult in ((32, 300, 60), (64, 300, 40), (128, 200, 25),
                       (256, 120, 16), (512, 40, 10)):
    tmax = mult * N * N
    lt, lk, el = run(N, nruns, tmax, seed=1000 + N)
    med = float(np.median(lt[lk])) if lk.any() else float('nan')
    print(f"{N:>6} {nruns:>5} {tmax:>8} {lk.sum():>8} {lk.mean():>7.3f} "
          f"{med:>14.0f} {med/N**2:>11.4f} {el:>7.1f}")
    sys.stdout.flush()
