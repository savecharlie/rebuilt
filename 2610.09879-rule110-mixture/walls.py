#!/usr/bin/env python3
"""walls.py -- the asymptotic velocity of ONE isolated domain wall.

WHY NOT A BIT FLIP. A single flip creates a PAIR of walls together, and the pair
carries a long transient: measured over successive windows out to t = 6000, the
ratio v_L/v_R falls 1.450 -> 1.245 -> 1.169 -> 1.140 -> 1.120 and is still
moving. A velocity read off any one window is a number about that window. The
paper quotes v_R = 0.219 and v_L = 0.236 (a 7% gap); a flip-pair measured over
[0,500) gives 0.173 and 0.251 (a 45% gap) and over [4500,6000) gives 0.212 and
0.238 (12%). All four pairs are the same system measured at different ages.

So: BUILD THE WALL DIRECTLY. Splice a lattice whose left half has been evolved
to phase a and whose right half to phase b, far apart on a big ring. Each wall
is then isolated, in steady state, with no pair-creation history, and v is the
drift of a single phase boundary.

A wall's position is found by PHASE CONSISTENCY, not by illegality: each block
form is compatible with a set of phases (`100` with all three, `001` with none),
so the lattice's phase label is read off the unambiguous blocks and the wall is
where it changes. Illegality counts the ragged interior of the out-of-phase
domain too, which is why the flip measurement above needed a burn-in at all.
"""
import sys

import numpy as np

RULES = (36, 37, 5)
TAB = np.array([[(r >> k) & 1 for k in range(8)] for r in RULES], np.uint8)
ALPH = {0: {(0, 0, 0), (1, 0, 0)},
        1: {(0, 1, 0), (0, 1, 1), (1, 0, 0), (1, 0, 1)},
        2: {(1, 0, 0), (1, 1, 0), (1, 1, 1)}}
# COMPAT[code] = bitmask of phases this block form can belong to
COMPAT = np.zeros(8, np.uint8)
for ph, A in ALPH.items():
    for (a, b, c) in A:
        COMPAT[4 * a + 2 * b + c] |= (1 << ph)


def step(s):
    res = np.tile(np.array([0, 1, 2], np.int64), s.shape[-1] // 3)
    k = (np.roll(s, 1, -1).astype(np.uint8) << 2) | (s << 1) | np.roll(s, -1, -1)
    return TAB[res, k]


def encode(bits):
    o = np.zeros((len(bits), 3), np.uint8)
    o[np.asarray(bits) == 0, 0] = 1
    return o.reshape(-1)


def phase_label(s, t):
    """Per-block phase label in {0,1,2}, or -1 where ambiguous/illegal.

    A block at global time t belongs to phase q if its form is compatible with
    (t - q) % 3 ... but what we actually want is the LOCAL beat offset. Label
    the block by which of the three alphabets its form is EXCLUSIVELY in; blocks
    compatible with more than one (i.e. `100`) are -1 and carry no information.
    """
    b = s.reshape(-1, 3).astype(np.int64)
    c = COMPAT[4 * b[:, 0] + 2 * b[:, 1] + b[:, 2]]
    lab = np.full(len(c), -1, np.int8)
    for ph in range(3):
        lab[c == (1 << ph)] = (ph - t) % 3
    return lab


def splice(N, offset, rng):
    """Ring of N blocks: left half at beat 0, right half `offset` beats ahead.

    Returns the lattice at a moment when both halves are legal in their own
    phase, and the two wall positions (N//4 and 3N//4 by construction).
    """
    a = encode(rng.integers(0, 2, N).astype(np.uint8))
    b = encode(rng.integers(0, 2, N).astype(np.uint8))
    for _ in range(offset):
        b = step(b)
    s = a.copy()
    q0, q1 = 3 * (N // 4), 3 * (3 * N // 4)
    s[q0:q1] = b[q0:q1]
    return s


def track(s, T, off, w=31, jump=60.0):
    """Follow the two phase boundaries by a SMOOTHED MAJORITY, not a count.

    TWO EARLIER VERSIONS WERE WRONG AND BOTH PRODUCED PLAUSIBLE NUMBERS.
      v1 located each wall as "the label change nearest its construction point",
      which silently re-anchors onto a spurious change once the wall has moved;
      a wall IS a disordered region and emits a cluster of them. Symptom: the
      left-movers came out stable to 1% while the right-movers disagreed between
      reps by a factor of two.
      v2 located each wall by COUNTING unambiguous blocks of each phase. But the
      block form `100` is legal in all three phases and so carries no label, and
      its ABUNDANCE DIFFERS BY PHASE -- phase 0's alphabet is {000,100}, so about
      half of a phase-0 domain is unlabelled, against about a quarter of a
      phase-1 domain. The count is therefore a phase-dependent fraction of the
      extent, and v2 reported velocities that changed sign between windows.

    What settled it was rendering the label field and LOOKING at it
    (`labels_off1.png`, `labels_off2.png`): the boundaries are dead straight,
    and the only hard part is that ~41% of blocks are unlabelled noise between
    them. So: +1 where the block is the spliced phase, -1 where it is the
    reference phase, 0 where unlabelled; box-smooth over `w` blocks; the wall is
    where the smoothed signal crosses zero. Unlabelled blocks dilute the signal
    but cannot move the crossing.
    """
    N = len(s) // 3
    pa = np.empty(T)
    pb = np.empty(T)
    ker = np.ones(w) / w
    prev = [N / 4.0, 3.0 * N / 4.0]
    for t in range(T):
        lab = phase_label(s, t)
        ind = np.where(lab == off % 3, 1.0, np.where(lab == 0, -1.0, 0.0))
        # v4, and the three earlier versions are kept in the record because each
        # produced plausible numbers while being wrong:
        #   v1  nearest label-change to the CONSTRUCTION point -- re-anchors onto
        #       a spurious change once the wall has moved away.
        #   v2  counts of unambiguously-labelled blocks -- `100` is legal in all
        #       three phases and its abundance DIFFERS BY PHASE (phase 0's
        #       alphabet is {000,100}, so half a phase-0 domain is unlabelled,
        #       against a quarter of a phase-1 domain), so the count is a
        #       phase-dependent fraction of the extent. Velocities changed sign.
        #   v3  argmin/argmax of cumsum(ind - mean(ind)). Parameter-free and
        #       right most of the time -- and in 3 of 10 reps at N = 6000 it
        #       reported the RIGHT-hand wall moving left at -0.2405, i.e. it
        #       returned the other boundary, because the mean it subtracts can
        #       cross zero (the labelled fraction depends on the rule-110
        #       content, which drifts) and the cumsum then flips sign globally.
        #       A 3-in-10 silent role swap is exactly the kind of thing that
        #       averages into a believable wrong number: it turned a clean
        #       +0.199 into +0.069 +- 0.068.
        # v4 tracks by CONTINUITY from the construction, capping the per-step
        # jump. The walls move <= 0.25 blocks/step and the observed per-step
        # excursion never exceeds 30, so a cap of `jump` cannot skip a wall.
        sm = np.convolve(np.r_[ind[-w:], ind, ind[:w]], ker, "same")[w:w + N]
        cr = np.flatnonzero(sm[:-1] * sm[1:] < 0).astype(float)
        for i in (0, 1):
            if len(cr):
                d = np.abs(cr - prev[i])
                d = np.minimum(d, N - d)            # the ring
                j = int(np.argmin(d))
                if d[j] <= jump:
                    prev[i] = float(cr[j])
            (pa if i == 0 else pb)[t] = prev[i]
        s = step(s)
    return pa, pb


def gate(rng):
    """A spliced lattice must be legal in its own phase away from the walls."""
    N = 600
    for off in (1, 2):
        s = splice(N, off, rng)
        lab = phase_label(s, 0)
        interior_l = lab[20:N // 4 - 20]
        interior_r = lab[N // 4 + 20:3 * N // 4 - 20]
        # the label of a region evolved `off` steps from an encoded state is
        # off % 3 -- "how many beats AHEAD of the reference this region is".
        # The first version of this gate expected (-off)%3 and the gate caught
        # it; that is the whole reason it runs before anything prints.
        for name, seg, want in (("left", interior_l, 0),
                                ("right", interior_r, off % 3)):
            v = seg[seg >= 0]
            if len(v) < 5 or not (v == want).all():
                return (f"gate: {name} interior of an offset-{off} splice is not "
                        f"a single phase ({set(v.tolist())} vs {want})")
    return None


if __name__ == "__main__":
    rng = np.random.default_rng(20261014)
    why = gate(rng)
    assert why is None, why
    print("gate passes: a spliced lattice is two clean single-phase domains\n")
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 6000
    T = int(sys.argv[2]) if len(sys.argv) > 2 else 4000
    reps = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    print(f"N = {N} blocks, T = {T} steps, {reps} reps per offset")
    print("offset 1 = spliced domain is one beat AHEAD; offset 2 = one beat BEHIND\n")
    out = {}
    for off in (1, 2):
        for r in range(reps):
            s = splice(N, off, rng)
            pa, pb = track(s, T, off)
            for name, P in (("left-edge", pa), ("right-edge", pb)):
                m = ~np.isnan(P)
                # only the stretch before the wedge closes / wraps
                good = np.flatnonzero(m)
                if len(good) < T // 3:
                    continue
                lo, hi = good[0], good[-1]
                span = hi - lo
                use = slice(lo + span // 10, lo + span - span // 10)
                tt = np.arange(T)[use][~np.isnan(P[use])]
                xx = P[use][~np.isnan(P[use])]
                if len(tt) < 100:
                    continue
                v = np.polyfit(tt, xx, 1)[0]
                out.setdefault((off, name), []).append(v)
    print(f"{'offset':>6} {'edge':>11} {'reps':>5} {'v (blocks/step)':>18}  reading L->R")
    for (off, name), vs in sorted(out.items()):
        vs = np.array(vs)
        # which phase order this edge sits at
        order = ("behind->ahead" if (off == 1 and name == "left-edge") or
                 (off == 2 and name == "right-edge") else "ahead->behind")
        print(f"{off:>6} {name:>11} {len(vs):>5}  {vs.mean():+.5f} +- "
              f"{vs.std(ddof=1)/len(vs)**0.5:.5f}   {order}")
    L = np.concatenate([out[k] for k in out if
                        (k[0] == 1 and k[1] == "left-edge") or
                        (k[0] == 2 and k[1] == "right-edge")])
    R = np.concatenate([out[k] for k in out if
                        (k[0] == 1 and k[1] == "right-edge") or
                        (k[0] == 2 and k[1] == "left-edge")])
    print(f"\n  behind->ahead walls (all move LEFT):  v_L = {abs(L.mean()):.5f} "
          f"+- {L.std(ddof=1)/len(L)**0.5:.5f}   n={len(L)}")
    print(f"  ahead->behind walls (all move RIGHT): v_R = {abs(R.mean()):.5f} "
          f"+- {R.std(ddof=1)/len(R)**0.5:.5f}   n={len(R)}")
    print(f"  paper: v_L 0.236, v_R 0.219 (+7.8%).  measured gap: "
          f"{100*(abs(L.mean())/abs(R.mean())-1):+.1f}%")
