#!/usr/bin/env python3
"""budgap.py -- P23: does a ring's gap drift track its OWN bud rate?

REGISTERED BEFORE RUNNING (16 Oct 2026, fire 344). Section 6 proposed that the
restoring force measured in section 4b -- a (2,2,2) triple relaxing toward equal
spacing and holding it, a (1,1,1) triple's smallest gap collapsing from 60 to 9 --
is not a separate force at all. It is what repeated budding-and-annihilating looks
like measured at the level of gaps: each escaped bud product moves a wall.

P23, with its sign fixed in advance:
  **Across rings, bud rate and gap drift are POSITIVELY correlated, and the
  correlation survives within the (1,1,1) chirality alone.**
The within-chirality clause is the part that matters. Across both chiralities a
positive correlation is nearly guaranteed, because (1,1,1) has both a thicker
cloud and a bigger drift -- that would be two group means masquerading as a
relationship. The test is Spearman rho within (1,1,1) only, n = 14 rings.

KILLED BY: rho <= 0 within (1,1,1), or a within-chirality p > 0.05 in a direction
I then explain away. If the across-group correlation is strong and the
within-group one is absent, the honest reading is that chirality sets both and
budding does not drive drift -- and I will say exactly that.

Gaps are measured from the THREE LONGEST plateaus, which are the persistent walls;
the buds are by construction the short ones. Drift is the standard deviation of
the three gaps, as section 4b used (it reported spread 90 -> 46 and 91 -> 240).
"""
import json

import numpy as np

from split import plateaus
from tri import Kit
from wallcount import locked_ring as wl
from wallcount import step as ws

kit = Kit((36, 37, 5))
N, T, EVERY = 300, 40_000, 25


def ring(beats, seed, gaps=(60, 90, 150)):
    """Section 4b's construction: three walls at deliberately unequal gaps."""
    rng = np.random.default_rng(seed)
    s = wl(N, rng)
    pos = np.cumsum([0] + list(gaps))[:3]
    for i, off in enumerate(beats):
        b = wl(N, rng)
        for _ in range(off):
            b = ws(b)
        lo = int(pos[i])
        hi = int(pos[i + 1]) if i + 1 < 3 else N
        s[3 * lo:3 * hi] = b[3 * lo:3 * hi]
    return s


def run(beats, seed):
    s = ring(beats, seed)
    buds, frames, spreads, alive = 0, 0, [], True
    for t in range(T):
        if t % EVERY == 0:
            pl = plateaus(kit, s, t)
            if len(pl) < 2:
                alive = False
                break
            frames += 1
            if len(pl) > 3:
                buds += 1
            top = sorted(pl, key=lambda r: -r[1])[:3]
            if len(top) == 3:
                st = sorted(r[2] for r in top)
                g = [st[1] - st[0], st[2] - st[1], N - (st[2] - st[0])]
                spreads.append((t, float(np.std(g))))
        s = kit.step(s)
    if not alive or len(spreads) < 50:
        return None
    sp = np.array(spreads)
    early = sp[sp[:, 0] < 4_000, 1].mean()
    late = sp[sp[:, 0] > T - 8_000, 1].mean()
    return dict(bud_rate=buds / frames, spread_early=float(early),
                spread_late=float(late), drift=float(late - early),
                frames=frames)


def spearman(x, y):
    rx = np.argsort(np.argsort(x)).astype(float)
    ry = np.argsort(np.argsort(y)).astype(float)
    rx -= rx.mean()
    ry -= ry.mean()
    return float((rx @ ry) / np.sqrt((rx @ rx) * (ry @ ry)))


if __name__ == "__main__":
    print("budgap.py -- P23, sign registered before running: rho > 0 WITHIN (1,1,1)")
    print(f"N={N}, gaps built (60,90,150), t to {T:,}, sampled every {EVERY}\n")
    out = {}
    for beats, lab in (((0, 1, 2), "(1,1,1)"), ((0, 2, 1), "(2,2,2)")):
        rows = [(sd, run(beats, 9000 + sd)) for sd in range(14)]
        good = [(sd, r) for sd, r in rows if r]
        print(f"=== {lab}: {len(good)}/14 rings still had three walls at t={T:,} ===")
        print(f"   {'seed':>4} {'bud rate':>9} {'spread 0-4k':>12} "
              f"{'spread late':>12} {'drift':>8}")
        for sd, r in good:
            print(f"   {sd:>4} {r['bud_rate']:>9.3f} {r['spread_early']:>12.1f} "
                  f"{r['spread_late']:>12.1f} {r['drift']:>+8.1f}")
        if len(good) >= 5:
            b = np.array([r['bud_rate'] for _, r in good])
            d = np.array([r['drift'] for _, r in good])
            rho = spearman(b, d)
            print(f"   Spearman rho(bud rate, drift) = {rho:+.3f}   n={len(good)}")
        out[lab] = [dict(seed=sd, **r) for sd, r in good]
    both = [(r['bud_rate'], r['drift']) for v in out.values() for r in v]
    if len(both) >= 6:
        b = np.array([x for x, _ in both])
        d = np.array([y for _, y in both])
        print(f"\nacross BOTH chiralities: rho = {spearman(b, d):+.3f}  n={len(both)}")
        print("(the across-group number is the one that can be two group means"
              " wearing a correlation's clothes -- the within-group one is the test)")
    with open("budgap.json", "w") as f:
        json.dump(out, f, indent=1)
