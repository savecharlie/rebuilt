#!/usr/bin/env python3
"""freevisit.py -- are the 5 never-visited table entries visited AT A WALL?

P14, registered in SECTION5_REGISTERED.md.

Exactly 19 of the 24 (residue, neighbourhood) entries are forced by the
emulation claim; the other 5 are never reached from any legal codeword
configuration (reachable.py, exhaustive to saturation).  A domain wall is not a
legal codeword configuration -- each side is, the junction is not.  If the
junction reaches a free entry, then the 32 main-cube triples, which agree on
every forced entry and may disagree on the free ones, are ONE system on the
computation and up to 32 different systems at the defects.

GATE (validated against a known answer before the new measurement speaks):
a single locked domain with no wall must visit ZERO free entries, because
reachable.py proves it exhaustively.  If the gate fires, the counter is wrong.
"""
import itertools
import sys

import numpy as np

from wallcount import step, walls, labels, locked_ring

RES = (0, 1, 2)
BASE = (36, 37, 5)


def _table(r):
    return {(k >> 2 & 1, k >> 1 & 1, k & 1): (r >> k) & 1 for k in range(8)}


def _nbr(s, i):
    n = len(s)
    return (s[(i - 1) % n], s[i], s[(i + 1) % n])


def _step_t(s, tabs):
    return tuple(tabs[i % 3][_nbr(s, i)] for i in range(len(s)))


def _enc(bits):
    out = []
    for b in bits:
        out += [0, 0, 0] if b else [1, 0, 0]
    return tuple(out)


def free_entries():
    """The 5 (residue, neighbourhood) pairs no legal configuration reaches."""
    tabs = {r: _table(BASE[r]) for r in RES}
    prev, vis = None, None
    for L in range(2, 11):
        vis = set()
        for bits in itertools.product((0, 1), repeat=L):
            s = _enc(bits)
            for _ in range(3):
                for i in range(len(s)):
                    vis.add((i % 3, _nbr(s, i)))
                s = _step_t(s, tabs)
        if vis == prev:
            break
        prev = vis
    allp = {(r, k) for r in RES for k in itertools.product((0, 1), repeat=3)}
    return sorted(allp - prev), sorted(prev)


FREE, FORCED = free_entries()
FREEMASK = np.zeros((3, 8), bool)
for r, k in FREE:
    FREEMASK[r, 4 * k[0] + 2 * k[1] + k[2]] = True


def visited_codes(s):
    """-> (n,) int array of 4*left+2*center+right for every site, and residues."""
    k = (np.roll(s, 1, -1).astype(np.int64) << 2) | (s.astype(np.int64) << 1) \
        | np.roll(s, -1, -1).astype(np.int64)
    res = np.tile(np.array([0, 1, 2]), s.shape[-1] // 3)
    return res, k


def run(charge, reps=10, T=600, burn=0, N=1800, seed=343, Rwin=40):
    """Count free-entry visits, and their block-offset from the nearest wall."""
    rng = np.random.default_rng(seed)
    hits = np.zeros((3, 8), np.int64)
    sites = 0
    off_hist = np.zeros(2 * Rwin + 1, np.int64)
    off_tot = 0
    for _ in range(reps):
        a = locked_ring(N, rng)
        b = locked_ring(N, rng)
        for _ in range(charge):
            b = step(b)
        s = a.copy()
        if charge:
            s[3 * (N // 4):3 * (3 * N // 4)] = b[3 * (N // 4):3 * (3 * N // 4)]
        for t in range(T):
            res, k = visited_codes(s)
            hit = FREEMASK[res, k]
            hits[res[hit], k[hit]] += 1
            sites += s.shape[-1]
            if t >= burn and hit.any():
                ch, _ = walls(s, t, 15)
                p = np.flatnonzero(ch)
                if len(p) == (2 if charge else 0) and charge:
                    hb = np.flatnonzero(hit) // 3           # block index of each hit
                    d = (hb[:, None] - p[None, :] + N // 2) % N - N // 2
                    near = d[np.arange(len(hb)), np.abs(d).argmin(1)]
                    inw = np.abs(near) <= Rwin
                    off_hist[near[inw] + Rwin] += 1
                    off_tot += len(near)
            s = step(s)
    return hits, sites, off_hist, off_tot


if __name__ == "__main__":
    print(f"free entries (never visited by any legal codeword config): {len(FREE)}")
    for r, k in FREE:
        print(f"   residue {r} (rule {BASE[r]}), neighbourhood {''.join(map(str,k))}"
              f"   -> rules differing here: "
              f"{[x for x in ({0:(36,38,52,54),1:(37,39),2:(5,13,21,29)}[r]) ]}")

    print("\nGATE -- a single locked domain, no wall, must visit ZERO free entries")
    h0, n0, _, _ = run(charge=0, reps=10, T=600)
    print(f"   {h0.sum()} free-entry visits in {n0:,} site-steps"
          f"   -> {'PASS' if h0.sum() == 0 else 'FAIL'}")
    if h0.sum():
        sys.exit("gate failed: the counter sees free entries where reachable.py proves none")

    for charge in (1, 2):
        h, n, off, tot = run(charge=charge, reps=10, T=600)
        print(f"\n=== charge-{charge} wall, 10 seeds x 600 steps, ring 1800 blocks ===")
        print(f"   free-entry visits: {h.sum():,} in {n:,} site-steps "
              f"({h.sum()/n:.3e} per site-step)")
        for r, k in FREE:
            c = h[r, 4 * k[0] + 2 * k[1] + k[2]]
            print(f"     residue {r} nbhd {''.join(map(str,k))}: {c:>9,}"
                  f"   {'VISITED' if c else 'never'}")
        if tot:
            w = np.flatnonzero(off) - 40
            print(f"   of visits assignable to a wall: {off.sum():,}/{tot:,} within"
                  f" 40 blocks; offsets {w.min():+d}..{w.max():+d}, "
                  f"peak at {int(np.argmax(off))-40:+d}")
