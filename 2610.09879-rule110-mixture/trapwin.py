#!/usr/bin/env python3
"""trapwin.py -- is "never 4 or 5 walls" physics, or is it the window?

NOTES.md lines 370-371 say of the trapped state: "Survivors are 3 or 6, always,
never 4 or 5, because a same-charge set needs n = 0 (mod 3)." That claim is in
the public rebuild. It was measured by builttrap.py with **W = 15**, and
crosscheck.py's CHECK B found that on live multi-wall rings the w=15 count reads
3 on frames where the raw-label plateau count reads 5 or 7.

So the claim has to be re-measured on the EXACT configuration it was made about:
builttrap's three arcs at beats (0,1,2), which is three walls of charge +1.

Note what is and is not at stake. The Z3 argument constrains SAME-CHARGE sets,
and 1+1+2+2 is a MIXED set, which Z3 never forbade. So the theory is not in
danger. The empirical sentence "never 4 or 5" is.

GATE 0, first, because two dynamics are being compared: wallcount.step and
tri.Kit.step must produce identical trajectories, or nothing below compares
anything.
"""
import numpy as np

from split import plateaus
from tri import Kit
from wallcount import locked_ring as wc_locked
from wallcount import step as wc_step
from wallcount import walls as wc_walls

kit = Kit((36, 37, 5))
N, W = 900, 15


def ring_of_beats(beats, seed):
    rng = np.random.default_rng(seed)
    s = wc_locked(N, rng)
    k = len(beats)
    for i, off in enumerate(beats):
        b = wc_locked(N, rng)
        for _ in range(off):
            b = wc_step(b)
        lo, hi = i * N // k, (i + 1) * N // k
        s[3 * lo:3 * hi] = b[3 * lo:3 * hi]
    return s


print("GATE 0 -- wallcount.step and tri.Kit.step must be the same dynamics")
s = wc_locked(N, np.random.default_rng(1))
a, b = s.copy(), s.copy()
same = True
for _ in range(300):
    a, b = wc_step(a), kit.step(b)
    if not np.array_equal(a, b):
        same = False
        break
print(f"   300 steps, identical: {same}  {'PASS' if same else 'FAIL'}\n")
assert same

print("builttrap's own configurations, counted BOTH ways, 8 seeds, t to 20,000")
print("   beats      seed   w=15 count histogram        raw-plateau histogram")
for beats, label in (((0, 1, 2), "(1,1,1)"), ((0, 2, 1), "(2,2,2)")):
    for sd in range(8):
        s = ring_of_beats(beats, 500 + sd)
        h15, hpl = {}, {}
        for t in range(20_000):
            if t >= 2_000 and t % 25 == 0:
                n15 = int(wc_walls(s, t, W)[1])
                npl = len(plateaus(kit, s, t))
                npl = 0 if npl < 2 else npl
                h15[n15] = h15.get(n15, 0) + 1
                hpl[npl] = hpl.get(npl, 0) + 1
            s = kit.step(s)
        f = lambda h: " ".join(f"{k}:{100*v/sum(h.values()):.0f}%"
                               for k, v in sorted(h.items()))
        print(f"   {beats} {label} {sd}   {f(h15):<26} {f(hpl)}")
    print()
