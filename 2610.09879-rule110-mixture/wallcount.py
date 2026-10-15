#!/usr/bin/env python3
"""wallcount.py -- how many domain walls are on the ring, and what they do.

Section 3 of arXiv:2610.09879. Predictions registered in SECTION3_REGISTERED.md
BEFORE this file existed.

A wall is a jump in the local BEAT PHASE. Phases live in Z3 and the total jump
around a periodic ring is 0 mod 3, so a wall carries a CHARGE in {+1,+2}, not a
sign. Two same-charge walls that meet give +2 == -1: they MERGE, 2 -> 1. Two
opposite-charge walls annihilate, 2 -> 0. That merge channel is why a generic
ring can reach zero walls at all; it is not in the paper.

COUNTING. phase_label() (from walls.py, the version that survived four bad
rulers) labels a block by the beat whose alphabet EXCLUSIVELY contains its form,
and -1 otherwise. About 41% of blocks are -1 -- the form `100` is legal in all
three beats -- so the label field is a sparse, noisy sampling of a piecewise
constant function. Fill it by voting: smooth the three phase indicators over a
box of w blocks and take the argmax. Count the places where the filled field
changes.

w IS A PARAMETER AND IT SETS THE ANSWER IN THE NOISY REGIME. G1 below pins it
against configurations whose wall number I built, across a range of w. Report
the span, never one value. (Three of my four bad rulers in fire 340 would have
been caught by a gate of this shape, and were not, because I did not build one.)
"""
import sys
import numpy as np

RULES = (36, 37, 5)
TAB = np.array([[(r >> k) & 1 for k in range(8)] for r in RULES], np.uint8)
ALPH = {0: {(0, 0, 0), (1, 0, 0)},
        1: {(0, 1, 0), (0, 1, 1), (1, 0, 0), (1, 0, 1)},
        2: {(1, 0, 0), (1, 1, 0), (1, 1, 1)}}
COMPAT = np.zeros(8, np.uint8)
for _ph, _A in ALPH.items():
    for (_a, _b, _c) in _A:
        COMPAT[4 * _a + 2 * _b + _c] |= (1 << _ph)
# EXCL[code] = the unique phase whose alphabet contains this form, else -1
EXCL = np.full(8, -1, np.int8)
for _code in range(8):
    m = COMPAT[_code]
    if m and (m & (m - 1)) == 0:
        EXCL[_code] = int(m).bit_length() - 1


def step(s):
    """One mixture step. s is (..., 3N) uint8."""
    res = np.tile(np.array([0, 1, 2], np.int64), s.shape[-1] // 3)
    k = (np.roll(s, 1, -1).astype(np.uint8) << 2) | (s << 1) | np.roll(s, -1, -1)
    return TAB[res, k]


def labels(s, t):
    """(..., 3N) uint8 -> (..., N) int8 beat label in {0,1,2}, -1 where ambiguous."""
    b = s.reshape(s.shape[:-1] + (s.shape[-1] // 3, 3)).astype(np.int64)
    code = 4 * b[..., 0] + 2 * b[..., 1] + b[..., 2]
    e = EXCL[code]
    out = np.where(e >= 0, (e - t) % 3, -1)
    return out.astype(np.int8)


def smooth_circ(x, w):
    """Circular box mean of odd width w along the last axis."""
    assert w % 2 == 1
    n = x.shape[-1]
    h = w // 2
    xx = np.concatenate([x[..., n - h:], x, x[..., :h]], axis=-1)
    c = np.cumsum(xx, axis=-1, dtype=np.float64)
    z = np.zeros(x.shape[:-1] + (1,))
    c = np.concatenate([z, c], axis=-1)
    return (c[..., w:] - c[..., :-w]) / w


def fill_nearest(lab):
    """Give every block the label of the nearest LABELLED block (circular).

    WHY THIS EXISTS. The first version smoothed the sparse label field directly
    and took the argmax. About 41% of blocks carry no label (the form `100` is
    legal in all three beats), so a chance run of w unlabelled blocks makes all
    three box sums zero and the argmax picks phase 0 arbitrarily -- a SPURIOUS
    WALL PAIR. Measured in resolution.py: the false-positive rate grew linearly
    with domain width and fell like 2**-w, so w=9 was exact at D=16 and wrong at
    D=20. Filling first removes the holes, so the smoothing that follows never
    sees a tie. Ties are the whole failure.
    """
    N = lab.shape[-1]
    three = np.concatenate([lab, lab, lab], axis=-1)
    idx = np.arange(3 * N)
    valid = three >= 0
    fwd = np.maximum.accumulate(np.where(valid, idx, -3 * N), axis=-1)
    bwd = np.minimum.accumulate(np.where(valid, idx, 6 * N)[..., ::-1], axis=-1)[..., ::-1]
    use_f = (idx - fwd) <= (bwd - idx)
    pick = np.where(use_f, fwd, bwd)
    pick = np.clip(pick, 0, 3 * N - 1)
    out = np.take_along_axis(three, pick, axis=-1)
    return out[..., N:2 * N]


def phase_field(s, t, w):
    """Filled, mode-filtered beat-phase field, (..., N) int8 in {0,1,2}.

    Fill the holes by nearest label, then take the majority over a box of w
    blocks. w is a RESOLUTION, not a tuning knob: domains narrower than about w
    are erased, and that limit is measured in resolution.py.
    """
    f = fill_nearest(labels(s, t).astype(np.int64))
    ind = np.stack([(f == p).astype(np.float64) for p in range(3)], axis=0)
    return smooth_circ(ind, w).argmax(0).astype(np.int8)


def walls(s, t, w):
    """-> (charge field, count).  charge[j] = (phase[j+1]-phase[j]) % 3, 0 = no wall."""
    ph = phase_field(s, t, w).astype(np.int16)
    ch = (np.roll(ph, -1, -1) - ph) % 3
    return ch, (ch != 0).sum(-1)


def encode(bits):
    o = np.zeros((len(bits), 3), np.uint8)
    o[np.asarray(bits) == 0, 0] = 1
    return o.reshape(-1)


def locked_ring(N, rng):
    return encode(rng.integers(0, 2, N).astype(np.uint8))


def spliced(N, arcs, offsets, rng):
    """Ring of N blocks with `arcs` out-of-phase arcs, each ahead by offsets[i] beats.
    Arcs are placed evenly and occupy a quarter of their slot. 2*arcs walls."""
    s = locked_ring(N, rng)
    for i in range(arcs):
        b = locked_ring(N, rng)
        for _ in range(offsets[i]):
            b = step(b)
        lo = (2 * i + 1) * N // (2 * arcs + 1)
        hi = lo + N // (2 * arcs + 2)
        s[3 * lo:3 * hi] = b[3 * lo:3 * hi]
    return s


def gate1(verbose=True):
    """G1: the counter returns the number of walls I built.

    THE GATE ENCODES THE MEASURED RESOLUTION LAW, not a hope. resolution.py
    (600-block ring, 40 random rings per cell, both offsets) shows the counter is
    exact for every domain width D >= w and erases domains narrower than w -- a
    clean monotone limit with no upper edge, once the label holes are filled
    before smoothing. So the gate requires EXACT for w <= D, and only REPORTS
    what happens for w > D. A gate that demanded exactness beyond the resolution
    would be a gate I had to tune away from, which is how bad rulers survive.
    """
    rng = np.random.default_rng(20261014)
    ws = [7, 11, 15, 21, 31, 41]
    cases = []
    for N in (128, 512):
        cases.append((f"locked N={N}", locked_ring(N, rng), 0, N))
        for off in (1, 2):
            cases.append((f"1 arc off={off} N={N}", spliced(N, 1, [off], rng), 2, N // 4))
        cases.append((f"2 arcs off=1,2 N={N}", spliced(N, 2, [1, 2], rng), 4, N // 6))
        cases.append((f"3 arcs off=1,2,1 N={N}", spliced(N, 3, [1, 2, 1], rng), 6, N // 8))
    ok = True
    for name, s, want, D in cases:
        got = [int(walls(s, 0, w)[1]) for w in ws]
        inres = [g for g, w in zip(got, ws) if w <= D]
        good = all(g == want for g in inres)
        ok &= good
        mark = "".join("." if w <= D else "~" for w in ws)
        if verbose:
            print(f"  {'ok ' if good else 'FAIL'} {name:<26} D={D:>4} want {want}  "
                  f"got {got}  in-res[{mark}]")
    s = spliced(256, 4, [1, 2, 1, 2], rng)
    got = int(walls(s, 0, 15)[1])
    if verbose:
        print(f"  {'ok ' if got == 8 else 'FAIL'} negative control: 4 arcs, w=15   want 8  got {got}")
    ok &= (got == 8)
    if verbose:
        print("  ('~' = w exceeds the domain width D, outside the instrument's "
              "resolution; not required to be exact)")
    return ok


if __name__ == "__main__":
    print("G1 -- counter against configurations whose wall number I built:")
    good = gate1()
    print("G1", "PASSED" if good else "FAILED")
    sys.exit(0 if good else 1)
