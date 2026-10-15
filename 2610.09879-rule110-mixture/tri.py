#!/usr/bin/env python3
"""tri.py -- the wall-measuring kit, parameterised by WHICH of the 42 emulating
triples is running.

WHY THIS EXISTS. freevisit.py (P14) measured that all five never-visited table
entries ARE visited at a domain wall, and nowhere else: zero in 32.4M site-steps
of wall-free lattice, 2.9e-4 per site-step with a wall, every assignable visit
within a dozen blocks of a wall centre.  The free entries are the entries where
the 42 emulating triples disagree.  So the triple is a knob on the DEFECT that
does not touch the COMPUTATION, and every instrument in this directory hardcodes
(36, 37, 5).

Everything here is wallcount.py's logic with RULES lifted to an argument, plus
two gates that must pass before a non-baseline triple is allowed to report:

  E -- the triple emulates rule 110 exhaustively for L <= 9 (so it is one of the
       42 and the computation is untouched);
  A -- its locked-domain ALPHABETS are identical to the baseline's, so the phase
       labeller -- which is derived from those alphabets -- is the same
       instrument and the numbers are comparable.  If A ever fails, the
       comparison is between two rulers and not between two walls.
"""
import itertools

import numpy as np

BASE = (36, 37, 5)
# the four families of 42, from NOTES.md section 3
MAIN_CUBE = [(a, b, c) for a in (36, 38, 52, 54) for b in (37, 39)
             for c in (5, 13, 21, 29)]
OTHERS = ([(70, b, c) for b in (41, 105) for c in (37, 165)]
          + [(84, b, c) for b in (37, 165) for c in (97, 105)]
          + [(a, 9, 86) for a in (82, 114)])
ALL42 = MAIN_CUBE + OTHERS


def _tab(rules):
    return np.array([[(r >> k) & 1 for k in range(8)] for r in rules], np.uint8)


def _d(r):
    return {(k >> 2 & 1, k >> 1 & 1, k & 1): (r >> k) & 1 for k in range(8)}


def _nbr(s, i):
    n = len(s)
    return (s[(i - 1) % n], s[i], s[(i + 1) % n])


def _enc(bits):
    out = []
    for b in bits:
        out += [0, 0, 0] if b else [1, 0, 0]
    return tuple(out)


def _dec(blk):
    out = []
    for j in range(0, len(blk), 3):
        b = blk[j:j + 3]
        if b == (1, 0, 0):
            out.append(0)
        elif b == (0, 0, 0):
            out.append(1)
        else:
            return None
    return tuple(out)


def emulates(rules, Lmax=9):
    """GATE E -- three mixture steps == one rule-110 step, exhaustively."""
    tabs = {r: _d(rules[r]) for r in (0, 1, 2)}
    t110 = _d(110)
    for L in range(1, Lmax + 1):
        for bits in itertools.product((0, 1), repeat=L):
            s = _enc(bits)
            for _ in range(3):
                s = tuple(tabs[i % 3][_nbr(s, i)] for i in range(len(s)))
            if _dec(s) != tuple(t110[_nbr(bits, i)] for i in range(len(bits))):
                return False
    return True


def alphabets(rules, L=8, reps=3):
    """The set of block forms a LOCKED domain shows at each beat, measured by
    running the triple on encoded configurations rather than asserted."""
    tabs = {r: _d(rules[r]) for r in (0, 1, 2)}
    A = {0: set(), 1: set(), 2: set()}
    for bits in itertools.product((0, 1), repeat=L):
        s = _enc(bits)
        for t in range(3 * reps):
            for j in range(0, len(s), 3):
                A[t % 3].add(s[j:j + 3])
            s = tuple(tabs[i % 3][_nbr(s, i)] for i in range(len(s)))
    return {k: frozenset(v) for k, v in A.items()}


BASE_ALPH = alphabets(BASE)


class Kit:
    """Stepper + wall instruments for one emulating triple."""

    def __init__(self, rules, check=True):
        self.rules = tuple(rules)
        self.TAB = _tab(self.rules)
        if check:
            assert emulates(self.rules), f"GATE E failed: {self.rules} does not emulate 110"
            a = alphabets(self.rules)
            assert a == BASE_ALPH, f"GATE A failed: {self.rules} has different locked alphabets"
        self.COMPAT = np.zeros(8, np.uint8)
        for ph, A in BASE_ALPH.items():
            for (x, y, z) in A:
                self.COMPAT[4 * x + 2 * y + z] |= (1 << ph)
        self.EXCL = np.full(8, -1, np.int8)
        for c in range(8):
            m = int(self.COMPAT[c])
            if m and (m & (m - 1)) == 0:
                self.EXCL[c] = m.bit_length() - 1

    # ---- dynamics -------------------------------------------------------
    def step(self, s):
        res = np.tile(np.array([0, 1, 2], np.int64), s.shape[-1] // 3)
        k = (np.roll(s, 1, -1).astype(np.uint8) << 2) | (s << 1) | np.roll(s, -1, -1)
        return self.TAB[res, k]

    # ---- phase reading (identical machinery for every triple, by GATE A) --
    def labels(self, s, t):
        b = s.reshape(s.shape[:-1] + (s.shape[-1] // 3, 3)).astype(np.int64)
        code = 4 * b[..., 0] + 2 * b[..., 1] + b[..., 2]
        e = self.EXCL[code]
        return np.where(e >= 0, (e - t) % 3, -1).astype(np.int8)

    @staticmethod
    def _smooth(x, w):
        n = x.shape[-1]
        h = w // 2
        xx = np.concatenate([x[..., n - h:], x, x[..., :h]], axis=-1)
        c = np.cumsum(xx, axis=-1, dtype=np.float64)
        z = np.zeros(x.shape[:-1] + (1,))
        c = np.concatenate([z, c], axis=-1)
        return (c[..., w:] - c[..., :-w]) / w

    @staticmethod
    def _fill(lab):
        N = lab.shape[-1]
        three = np.concatenate([lab, lab, lab], axis=-1)
        idx = np.arange(3 * N)
        valid = three >= 0
        fwd = np.maximum.accumulate(np.where(valid, idx, -3 * N), axis=-1)
        bwd = np.minimum.accumulate(np.where(valid, idx, 6 * N)[..., ::-1], axis=-1)[..., ::-1]
        pick = np.clip(np.where((idx - fwd) <= (bwd - idx), fwd, bwd), 0, 3 * N - 1)
        return np.take_along_axis(three, pick, axis=-1)[..., N:2 * N]

    def phase_field(self, s, t, w):
        f = self._fill(self.labels(s, t).astype(np.int64))
        ind = np.stack([(f == p).astype(np.float64) for p in range(3)], axis=0)
        return self._smooth(ind, w).argmax(0).astype(np.int8)

    def walls(self, s, t, w):
        ph = self.phase_field(s, t, w).astype(np.int16)
        ch = (np.roll(ph, -1, -1) - ph) % 3
        return ch, (ch != 0).sum(-1)

    # ---- configurations -------------------------------------------------
    @staticmethod
    def encode(bits):
        o = np.zeros((len(bits), 3), np.uint8)
        o[np.asarray(bits) == 0, 0] = 1
        return o.reshape(-1)

    def locked_ring(self, N, rng):
        return self.encode(rng.integers(0, 2, N).astype(np.uint8))

    def two_domain(self, N, charge, rng):
        """Ring of N blocks: right half `charge` beats ahead. Exactly 2 walls."""
        a = self.locked_ring(N, rng)
        b = self.locked_ring(N, rng)
        for _ in range(charge):
            b = self.step(b)
        s = a.copy()
        s[3 * (N // 4):3 * (3 * N // 4)] = b[3 * (N // 4):3 * (3 * N // 4)]
        return s

    def ring3(self, N, charge, gaps, rng):
        """Three walls of equal charge at built positions sum(gaps[:k])."""
        pos = np.cumsum([0] + list(gaps[:-1]))
        s = self.locked_ring(N, rng)
        cur = 0
        for k in range(3):
            lo = int(pos[k])
            hi = int(pos[k + 1]) if k + 1 < 3 else N
            b = self.locked_ring(N, rng)
            cur = (cur + charge) % 3
            for _ in range(cur):
                b = self.step(b)
            s[3 * lo:3 * hi] = b[3 * lo:3 * hi]
        return s


if __name__ == "__main__":
    print(f"all 42 emulating triples, gates E and A:")
    bad = []
    for t in ALL42:
        try:
            Kit(t)
            ok = "ok"
        except AssertionError as e:
            ok = f"FAIL ({e})"
            bad.append(t)
        print(f"   {str(t):<18} {ok}")
    print(f"\n{len(ALL42) - len(bad)}/{len(ALL42)} pass both gates")
