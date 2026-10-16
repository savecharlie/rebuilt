#!/usr/bin/env python3
"""burnin.py -- the correction to split.py, and the thing it measures instead.

WHAT WENT WRONG, first, because it reached a results file.
`split.py::time_to_split` scanned for a third plateau from t = 0 on a ring made
by SPLICING two locked rings together. A splice is not a relaxed wall. It
reported median split times of 5 steps (charge 2) and 2 steps (charge 1) -- and
what it was timing was the splice relaxing into a wall, not a wall splitting.
Measuring from t=0 on a construction that is not a steady state measures the
CONSTRUCTION. That is the seventh instrument of mine to lie on this paper and
it is a new species: not a selecting filter and not a contaminating one, but a
TRANSIENT of the initial condition read as a dynamical event. Its signature is
an effect whose timescale equals the relaxation time and not any physical one.

WHAT IS MEASURED HERE, after a burn-in:
  1. the steady-state FREQUENCY of a third plateau at the intermediate phase,
     per charge -- which is what "does a wall split" actually means;
  2. the strip-length distribution when one is present;
  3. the same for the `arc` construction race.py used, to test fire 343's claim
     that split walls "sit five blocks apart and persist to t = 20,000".

The detector is split.py's, with the response function measured there:
P(detect) = 1.00 for a strip of >= 8 blocks at either charge, 0.78/0.53 at one
block, and a flat length bias of -1 (charge 2) / -2 (charge 1) blocks.
"""
import json
import sys

import numpy as np

from split import built_wall, plateaus, strip_len
from tri import Kit

BASE = (36, 37, 5)
BURN = 1_000          # steps discarded; the splice relaxes in O(10)
N = 600


def steady(kit, make, charge, seed, T=30_000, every=10):
    """Frequency and length distribution of a third plateau after burn-in."""
    host = {2: 1, 1: 2}[charge]
    s = make(kit, charge, seed)
    for _ in range(BURN):
        s = kit.step(s)
    seen, frames, lens, offphase, nplat = 0, 0, [], 0, []
    for t in range(BURN, BURN + T):
        if t % every == 0:
            pl = plateaus(kit, s, t)
            frames += 1
            nplat.append(len(pl))
            if len(pl) >= 3:
                sh = min(pl, key=lambda r: r[1])
                if sh[0] == host:
                    seen += 1
                    lens.append(sh[1])
                else:
                    offphase += 1
        s = kit.step(s)
    return dict(frames=frames, seen=seen, freq=seen / frames,
                offphase=offphase, lens=lens,
                mean_nplat=float(np.mean(nplat)),
                max_nplat=int(max(nplat)))


def mk_built(kit, charge, seed):
    return built_wall(kit, N, charge, 0, seed)


def mk_arc(kit, charge, seed, d=200):
    """race.py's construction: an arc of d blocks `charge` beats ahead."""
    rng = np.random.default_rng(seed)
    s = kit.locked_ring(N, rng)
    b = kit.locked_ring(N, np.random.default_rng(seed))
    for _ in range(charge):
        b = kit.step(b)
    lo = 200
    s[3 * lo:3 * (lo + d)] = b[3 * lo:3 * (lo + d)]
    return s


if __name__ == "__main__":
    kit = Kit(BASE)
    print("burnin.py -- steady-state frequency of a third plateau")
    print(f"triple {BASE}, N={N} blocks, burn-in {BURN:,} steps, "
          f"then 30,000 steps sampled every 10\n")
    out = {}
    for name, make in (("built (spliced wall)", mk_built),
                       ("arc d=200 (race.py's)", mk_arc)):
        print(f"=== {name} ===")
        for charge in (2, 1):
            fr, ln, off = [], [], 0
            for sd in range(6):
                r = steady(kit, make, charge, 300 + sd)
                fr.append(r['freq'])
                ln += r['lens']
                off += r['offphase']
            fr = np.array(fr)
            ln = np.array(ln)
            bias = {2: 1, 1: 2}[charge]
            print(f"  charge {charge}: third plateau in "
                  f"{100*fr.mean():.2f}% of frames "
                  f"(per-seed {', '.join(f'{100*x:.1f}' for x in fr)})")
            if len(ln):
                print(f"            when present: raw length mean {ln.mean():.2f} "
                      f"median {np.median(ln):.0f} max {ln.max()}  "
                      f"-> bias-corrected mean {ln.mean()+bias:.2f}")
                u, c = np.unique(ln, return_counts=True)
                print(f"            raw length histogram "
                      f"{dict(zip(u.tolist(), c.tolist()))}")
            print(f"            off-phase third plateaus: {off}")
            out[f"{name}|{charge}"] = dict(
                freq=float(fr.mean()), freq_per_seed=fr.tolist(),
                n=int(len(ln)),
                raw_mean=float(ln.mean()) if len(ln) else None,
                corrected_mean=float(ln.mean() + bias) if len(ln) else None,
                offphase=off)
        print()
    with open("burnin.json", "w") as f:
        json.dump(out, f, indent=1)
    print("wrote burnin.json")
