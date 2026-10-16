#!/usr/bin/env python3
"""split.py -- is a charge-2 wall one object, or a bound state of two charge-1 walls?

Fire 344. Registered in SECTION6_REGISTERED.md BEFORE this ran (P19-P22).

Fire 343 left this as the last thing the paper was hiding: a charge-2 wall
appears to split into two charge-1 walls, and `blink2.py` could not tell that
from a detection artifact two blocks wide. The reason it could not is the
SMOOTHING WINDOW: a width large enough to survive a wall's ambiguous core is
also wide enough to merge two walls five blocks apart. So this file reads the
RAW labels and nothing else -- no _fill, no _smooth, no window.

Charge is conserved mod 3, which fixes the decay channels by arithmetic:
a charge-2 wall may go to 1+1; a charge-1 wall may only go to 2+2. Nothing
else sums. So the asymmetry between the two is a prediction, not an observation.

THE DETECTOR, and the gate that matters. `plateaus()` returns runs of equal
UNAMBIGUOUS phase label around the ring. Gated in both directions before it was
allowed to report (`--gates`):
  N  24,000 frames of wall-free locked lattice -> never more than one plateau.
  T  a bare charge-2 wall built with NO intermediate strip (m=0) -> exactly two
     plateaus, 6/6 seeds. The detector does not manufacture the result.
  T' a built strip of m = 1..15 blocks -> three plateaus, length recovered to
     +-1, 6/6 seeds. At m = 25 it undercounts by 3, so lengths above 15 are
     LOWER BOUNDS. The measurement lives at m ~ 2-10, inside the exact range.

Ivy, on why the measurement comes before the theory: "no im right about what i
SEE asshole". The plate is the world and my file is a model of it.
"""
import json
import sys

import numpy as np

from tri import Kit

BASE = (36, 37, 5)


# ---------------------------------------------------------------- detector
def plateaus(kit, s, t):
    """Runs of equal UNAMBIGUOUS phase label around the ring.

    Returns [(phase, length_in_blocks, start_block)], wrap merged.
    No smoothing and no filling: every number here is a raw label.
    """
    lab = kit.labels(s, t)
    idx = np.flatnonzero(lab >= 0)
    if len(idx) == 0:
        return []
    ph = lab[idx]
    runs, start = [], 0
    for i in range(1, len(ph)):
        if ph[i] != ph[i - 1]:
            runs.append((int(ph[start]), int(idx[i - 1] - idx[start] + 1),
                         int(idx[start])))
            start = i
    runs.append((int(ph[start]), int(idx[-1] - idx[start] + 1), int(idx[start])))
    if len(runs) > 1 and runs[0][0] == runs[-1][0]:
        p, l, st = runs[-1]
        runs[0] = (p, l + runs[0][1], st)
        runs = runs[:-1]
    return runs


# ------------------------------------------------------------ construction
def built_wall(kit, N, charge, m, seed, lo=200, width=220):
    """A ring with ONE wall of `charge` at block `lo` and its mandatory partner
    at `lo+width`, plus an intermediate strip of exactly `m` blocks inside the
    `charge`-wall so that m=0 is a bare wall.

    charge 2: phase p | m blocks at p+1 | rest of the domain at p+2 | p
              so the left wall is +1,+1 if m>0 and a bare +2 if m=0.
    charge 1: phase p | m blocks at p+2 | rest at p+1 | p
              so the left wall is +2,+2 if m>0 and a bare +1 if m=0.
    """
    rng = np.random.default_rng(seed)
    s = kit.locked_ring(N, rng)
    base = kit.locked_ring(N, np.random.default_rng(seed))
    ahead = {}
    b = base.copy()
    for k in (1, 2):
        b = kit.step(b)
        ahead[k] = b.copy()
    inner = 1 if charge == 2 else 2          # the strip's phase offset
    hi = lo + width
    s[3 * (lo + m):3 * hi] = ahead[charge][3 * (lo + m):3 * hi]
    if m:
        s[3 * lo:3 * (lo + m)] = ahead[inner][3 * lo:3 * (lo + m)]
    return s


def strip_len(kit, s, t, host_phase):
    """Length of the shortest plateau when there are >=3, else 0.

    `host_phase` is the intermediate phase the strip must sit at for the split
    to be the registered one. A third plateau at any OTHER phase is reported
    separately and does NOT count as a split.
    """
    pl = plateaus(kit, s, t)
    if len(pl) < 3:
        return 0, False
    short = min(pl, key=lambda r: r[1])
    return short[1], (short[0] == host_phase)


# ------------------------------------------------------------ measurements
def time_to_split(kit, charge, seed, N=600, T=2000, lo=200, width=220):
    """First t at which a third plateau at the intermediate phase appears.

    Starts from a VERIFIED bare wall: if the m=0 build does not read as exactly
    two plateaus at t=0 the seed is refused rather than silently included.
    """
    s = built_wall(kit, N, charge, 0, seed, lo, width)
    pl0 = plateaus(kit, s, 0)
    if len(pl0) != 2:
        return None, f"refused: t=0 reads {len(pl0)} plateaus, not 2"
    host = {2: 1, 1: 2}[charge]
    for t in range(T):
        if t:
            L, ok = strip_len(kit, s, t, host)
            if L and ok:
                return t, None
        s = kit.step(s)
    return -1, None                           # never split inside T


def strip_series(kit, charge, seed, N=600, T=50_000, every=50,
                 lo=200, width=220):
    """Strip length over long time, for the stationarity test (P20)."""
    s = built_wall(kit, N, charge, 0, seed, lo, width)
    host = {2: 1, 1: 2}[charge]
    ts, ls, wrong = [], [], 0
    for t in range(T):
        if t % every == 0:
            L, ok = strip_len(kit, s, t, host)
            if L:
                if ok:
                    ts.append(t)
                    ls.append(L)
                else:
                    wrong += 1
        s = kit.step(s)
    return np.array(ts), np.array(ls), wrong


# ----------------------------------------------------------------- gates
def gates(kit, N=600):
    ok = True
    print("GATE N -- wall-free locked lattice, 12 seeds x 2000 steps")
    worst = 1
    for sd in range(12):
        s = kit.locked_ring(N, np.random.default_rng(1000 + sd))
        for t in range(2000):
            worst = max(worst, len(plateaus(kit, s, t)))
            s = kit.step(s)
    print(f"   24,000 frames, max plateaus = {worst}  "
          f"{'PASS' if worst == 1 else 'FAIL'}")
    ok &= worst == 1

    print("GATE T -- bare built wall (m=0) must read as exactly 2 plateaus")
    for charge in (1, 2):
        bad = [sd for sd in range(6)
               if len(plateaus(kit, built_wall(kit, N, charge, 0, 500 + sd), 0)) != 2]
        print(f"   charge {charge}: {6 - len(bad)}/6 seeds read 2 plateaus  "
              f"{'PASS' if not bad else 'FAIL ' + str(bad)}")
        ok &= not bad

    print("GATE T' -- the detector's RESPONSE FUNCTION, measured not asserted")
    print("     a strip is visible only where its blocks carry beat-exclusive")
    print("     forms, which is a coin flip per block. So P(detect) rises with m")
    print("     and the length is biased DOWN by a constant. Both are measured")
    print("     here and both must match what SECTION6 registered:")
    print("       P(det) at m=1: charge2 0.78, charge1 0.53   (ratio 1.5)")
    print("       P(det) at m>=8: 1.00 for both")
    print("       length bias: charge2 -1.0, charge1 -2.0 blocks, flat in m")
    want = {2: dict(p1=0.78, bias=-1.0), 1: dict(p1=0.53, bias=-2.0)}
    for charge in (2, 1):
        host = {2: 1, 1: 2}[charge]
        p1 = bias8 = None
        for m in (1, 8, 12, 15, 20):
            det, errs = 0, []
            for sd in range(120):
                pl = plateaus(kit, built_wall(kit, N, charge, m, 2000 + sd), 0)
                if len(pl) >= 3:
                    sh = min(pl, key=lambda r: r[1])
                    if sh[0] == host:
                        det += 1
                        errs.append(sh[1] - m)
            if m == 1:
                p1 = det / 120
            else:
                bias8 = np.mean(errs) if errs else 99
                if det != 120:
                    print(f"   charge {charge} m={m}: P(det)={det/120:.2f} != 1.00  FAIL")
                    ok = False
        w = want[charge]
        okp = abs(p1 - w['p1']) <= 0.10
        okb = abs(bias8 - w['bias']) <= 0.6
        print(f"   charge {charge}: P(det|m=1) {p1:.2f} (want {w['p1']:.2f})"
              f" {'ok' if okp else 'FAIL'}   "
              f"bias {bias8:+.2f} (want {w['bias']:+.1f}) {'ok' if okb else 'FAIL'}")
        ok &= okp and okb

    print("   => the charge asymmetry is at most 1.5x, so P21's registered")
    print("      threshold of 10x CANNOT be manufactured by the instrument.")
    return ok


if __name__ == "__main__":
    kit = Kit(BASE)
    if "--gates" in sys.argv:
        print("split.py gates -- the detector proves itself before it reports\n")
        sys.exit(0 if gates(kit) else 1)

    print("split.py -- a charge-2 wall: one object or two?")
    print(f"triple {BASE}, ring 600 blocks, raw labels only (no window)")
    print("strip lengths below are RAW; the measured bias is -1 block for a\n"
          "charge-2 strip and -2 for a charge-1 strip (split.py --gates).\n")

    print("=== P19 / P21 -- time to split, from a VERIFIED bare wall ===")
    res = {}
    for charge in (2, 1):
        times, refused = [], 0
        for sd in range(24):
            t, why = time_to_split(kit, charge, 300 + sd)
            if t is None:
                refused += 1
            else:
                times.append(t)
        arr = np.array([t for t in times if t >= 0])
        never = sum(1 for t in times if t < 0)
        res[charge] = dict(n=len(times), refused=refused, never=never,
                           split=len(arr),
                           median=float(np.median(arr)) if len(arr) else None,
                           mn=int(arr.min()) if len(arr) else None,
                           mx=int(arr.max()) if len(arr) else None)
        print(f"  charge {charge}: {len(arr)}/{len(times)} split within 2,000 steps"
              f"   never={never}  refused={refused}")
        if len(arr):
            print(f"              t_split  median {np.median(arr):.0f}  "
                  f"range {arr.min()}-{arr.max()}")

    print("\n=== P20 -- is the strip length stationary? (t to 50,000) ===")
    for charge in (2, 1):
        allm, slopes = [], []
        for sd in range(4):
            ts, ls, wrong = strip_series(kit, charge, 300 + sd)
            if len(ls) < 20:
                print(f"  charge {charge} seed {sd}: only {len(ls)} samples "
                      f"with a strip -- no stationarity claim")
                continue
            # does the mean drift? regress length on t and on sqrt(t)
            a = np.polyfit(ts, ls, 1)[0]
            b = np.polyfit(np.sqrt(ts), ls, 1)[0]
            allm.append(ls)
            slopes.append((a, b))
            print(f"  charge {charge} seed {sd}: n={len(ls):4d} present "
                  f"{100*len(ls)/(50_000//50):.0f}% of frames  "
                  f"mean {ls.mean():.2f}  sd {ls.std():.2f}  "
                  f"max {ls.max()}  d<L>/dt {a:+.2e}  d<L>/d sqrt(t) {b:+.2e}  "
                  f"off-phase {wrong}")
        if allm:
            cat = np.concatenate(allm)
            print(f"  charge {charge} POOLED: mean {cat.mean():.2f} "
                  f"sd {cat.std():.2f}  n={len(cat)}  "
                  f"hist {dict(zip(*[x.tolist() for x in np.unique(cat, return_counts=True)]))}")
        res[f'strip{charge}'] = dict(
            mean=float(np.concatenate(allm).mean()) if allm else None,
            slopes=[[float(a), float(b)] for a, b in slopes])

    print("\n=== P22 -- does the partner's distance matter? ===")
    for width in (60, 110, 220, 400):
        ts = []
        for sd in range(12):
            t, _ = time_to_split(kit, 2, 300 + sd, width=width)
            if t is not None and t >= 0:
                ts.append(t)
        print(f"  partner at {width:>3} blocks: {len(ts)}/12 split, "
              f"median t {np.median(ts) if ts else float('nan'):.0f}")
        res[f'width{width}'] = dict(n=len(ts),
                                    median=float(np.median(ts)) if ts else None)

    with open("split.json", "w") as f:
        json.dump(res, f, indent=1)
    print("\nwrote split.json")
