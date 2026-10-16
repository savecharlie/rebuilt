#!/usr/bin/env python3
"""crosscheck.py -- wallseq's plateau count against TWO independent instruments.

THE GATE I WROTE THAT CANNOT FAIL. wallseq.py refuses a frame whose wall charges
do not sum to 0 mod 3, and I wrote in its docstring that this is "a free CHECK on
every frame". It is a TAUTOLOGY. The charges are ((p[i+1]-p[i]) mod 3) around a
cycle, and sum_i (p[i+1]-p[i]) telescopes to exactly 0, so the sum of the
residues is identically 0 mod 3. It can never fail, and it did not: refused = 0
in 96,000 frames, which I read as the instrument being healthy. Ninth ruler on
this paper, and the first one whose defect is that it is VACUOUSLY TRUE. The tell
I want to keep: a check that has never once fired is not evidence of health until
I have shown it CAN fire.

The two real checks, both instrument-against-instrument:

  A  LOCKED agreement. wallseq reads 1 plateau (0 walls) in 30 of 32 runs. An
     independent criterion -- race.py's, which touches no labels at all -- is
     that every block is a codeword in {0,4} and none is illegal. If plateaus
     say "locked" and the raw codes say "not locked", the plateau reading is
     wrong.
  B  WALL-COUNT agreement where both instruments are valid. On the live runs,
     compare the plateau count with tri.Kit.walls()'s smoothed count at a window
     small enough to resolve the spacing. They disagree near the resolution
     limit by construction; they must agree when walls are far apart.
"""
import numpy as np

from invade import ILLEGAL, codes
from split import plateaus
from tri import Kit
from wallseq import arc, wall_charges

kit = Kit((36, 37, 5))
BURN, N = 1_000, 600


def raw_locked(s):
    """race.py's criterion. No labels, no plateaus, no window."""
    c = codes(s)
    return (not np.isin(c, ILLEGAL).any()) and bool(np.isin(c, [0, 4]).all())


print("CHECK A -- plateau 'locked' vs the raw-code criterion (no labels at all)")
print(f"   {'cell':<18} {'seed':>4} {'plateaus':>9} {'raw locked':>11} "
      f"{'illegal frac':>13}  agree")
dis = 0
for ahead in (1, 2):
    for d in (200, 300):
        for sd in range(8):
            s = arc(kit, d, ahead, 300 + sd)
            for _ in range(BURN + 5_000):
                s = kit.step(s)
            t = BURN + 5_000
            npl = len(plateaus(kit, s, t))
            rl = raw_locked(s)
            fi = float(np.isin(codes(s), ILLEGAL).mean())
            ok = (npl == 1) == rl
            dis += (not ok)
            if sd < 3 or not ok:
                print(f"   d={d} ahead={ahead:<6} {sd:>4} {npl:>9} "
                      f"{str(rl):>11} {fi:>13.4f}  {'ok' if ok else 'DISAGREE'}")
print(f"   disagreements: {dis}/32  "
      f"{'PASS -- the two instruments agree' if dis == 0 else 'FAIL'}")

print("\nCHECK B -- plateau wall count vs the smoothed walls() count, live runs")
print("   (the two live cells are d=200 ahead=1, seeds 2 and 5)")
for sd in (2, 5):
    s = arc(kit, 200, 1, 300 + sd)
    for _ in range(BURN):
        s = kit.step(s)
    rows = []
    for k in range(12):
        t = BURN + k * 500
        ch, ph = wall_charges(kit, s, t)
        npl = len(ch) if ch else 0
        sm = {w: int(kit.walls(s, t, w)[1]) for w in (3, 5, 9, 15)}
        rows.append((t, npl, ''.join(map(str, ph)), sm))
        for _ in range(500):
            s = kit.step(s)
    print(f"   seed {sd}:   t     plateaus(charges)   smoothed walls at w=3/5/9/15")
    for t, npl, ph, sm in rows:
        print(f"           {t:>6}   {npl:>2} walls  [{ph}]      "
              f"{sm[3]:>3} {sm[5]:>3} {sm[9]:>3} {sm[15]:>3}")
