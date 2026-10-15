#!/usr/bin/env python3
"""tracked.py -- walls that keep their NAMES, and no way to ask for a sorted gap.

WHY THIS MODULE EXISTS. In fire 341 four separate instruments broke the same
way inside one session, and the fourth broke after I had already found and
written up the first three:

  * `lastpair.py` discarded any pair whose wall count changed during its
    measurement window -- which is exactly what a pair about to close does. It
    threw away 164 of 195 cases and shrank the effect from 19.7x to 11.8x.
  * `diffuse.py` took the variance of three gaps that are constrained to sum to
    N, so two of its three "independent" samples were determined by the first.
  * the first repulsion test averaged over runs contributing different numbers
    of frames.
  * `forcelaw.py` tracked *the minimum* of three fluctuating gaps and reported
    its order-statistic drift as a force law.

Four bugs, one mistake: **a statistic over a set whose membership depends on the
thing being measured.** Sorted, minimum, surviving, readable -- every one of
those words is a filter correlated with the answer. And fire 342's `nullwalk.py`
measured the two signs, in a system with no force at all: sorting makes three
gaps look like they are SPREADING (90 -> 104.6), survivorship makes them look
like they are EQUALISING (90 -> 49). Either artifact alone is big enough to be
the whole of a result.

A resolution to be careful does not survive my own reset. So this module makes
the bad version unavailable rather than discouraged:

  1. Walls are NAMED at construction and followed by continuity. The matching is
     a bijection, so two names can never collapse onto one detection -- the old
     `prev[i]` loop in walls.py could do that silently and did.
  2. The first change in the wall count is a HARD STOP, recorded as `t_event`.
     There is no filter applied afterwards, because there is no data afterwards.
     A censored mean cannot be computed from something that was never recorded.
  3. There is no function here that returns a sorted gap, a minimum gap, or a
     mean over a subset. `Track.gaps` is indexed by name.
  4. Identity gaps must sum to N at every frame, checked every frame. A tracker
     that swapped two walls gives 0 or 2N and raises instead of reporting.

A BLINK IS NOT AN EVENT, AND THE FIRST VERSION OF THIS MODULE CONFLATED THEM.
Stopping on any change in wall count ended charge-(1,1,1) tracks after a median
of 650 steps, against fire 341's median resolution of 54,900 for the same rings
-- a factor of 80 between two of my own measurements, which can only be the
instrument. `flicker.py` recorded the raw count without stopping: about 2% of
frames show a spurious EXTRA wall, roughly one excursion per 250 steps, and in
the runs examined every single one returned to 3 within one sample. `blinkwhere.py`
then located them: 100% of extra detections sit within 9 blocks of a real wall,
median 2, i.e. one wall momentarily resolved as two. So:

  * MORE detections than names is a blink. Match the names to the nearest
    detections (the assignment is injective on detections, so extras are simply
    unclaimed) and carry on, counting it in `self.blinks`.
  * FEWER detections than names is a candidate annihilation, and it has to
    PERSIST for `persist` samples before it is believed. If the count comes
    back, the track resumes and the skipped samples are recorded in
    `self.holes` -- never silently interpolated.

`t_event` is the time of the FIRST shortfall of a persistent run, not the time
the persistence test completed, so no lookahead is needed and the number means
what it says.

Gates (run `python3 tracked.py`): the sum invariant; recovery of the isolated
wall velocities measured independently in fire 341 (charge 1 at -0.2414, charge
2 at +0.1932, 74 sigma apart, no overlap in 38 walls); the stop guard firing on
a ring built to annihilate, at a time its own closing speed predicts; and the
blink guard surviving a (1,1,1) ring for thousands of steps where the first
version died at 650.
"""
import numpy as np

from wallcount import step, walls, locked_ring

W = 15                 # phase-field smoothing width, as used throughout section 3
MAXV = 0.25            # measured ceiling on |wall velocity|, blocks per step
# JITTER is the DETECTOR's per-sample position scatter, which is a different
# thing from the wall's motion and is much bigger. jump.py measures it over
# 83,000 sampled frames: at every=1, where the wall moves 0.25 blocks, the
# detection moves a median of 0 and a maximum of 8, because the position is a
# majority vote over a 15-block box on a label field that is ~40% holes. The
# first version of this module capped the jump at MAXV*every*3+2 = 5.75 blocks
# at every=5 and was thrown off healthy rings by real 6-block jumps. Measured
# maxima: 8 at every=1, 10 at every=5, 16 at every=25.
JITTER = 12.0          # blocks of slack for detector scatter, from jump.py


class Lost(Exception):
    """Tracking failed -- as distinct from a wall physically disappearing."""


class Track:
    """Named walls on a ring, followed until the wall count changes.

    names    tuple of (slot, charge) in cyclic order, fixed at t0
    t        sampled times, length = frames
    pos      (frames, n) unwrapped positions, float
    gaps     (frames, n) gaps[k, i] = blocks from name i to name i+1 at t[k]
    t_event  the time the wall count first differed from n, or None
    """

    def __init__(self, s, N, t0=0, w=W, every=1, persist=4):
        self.N, self.w, self.every, self._t0 = N, w, every, t0
        self.persist = persist
        self.blinks, self.holes = 0, []
        self.unmatched = 0
        ch, _ = walls(s, t0, w)
        p = np.flatnonzero(ch)
        if len(p) < 2:
            raise Lost(f"only {len(p)} walls at t0; nothing to name")
        self.charges = tuple(int(c) for c in ch[p])
        self.n = len(p)
        self._x = p.astype(float)          # unwrapped, starts equal to the ring index
        self._last = p.astype(float) % N
        self._pos, self._t, self._g = [self._x.copy()], [t0], [self._gaps()]
        self.t_event = None
        # published immediately so the initial gaps can be read (and a ring
        # rejected) before any stepping happens
        self.t = np.array(self._t)
        self.pos = np.array(self._pos)
        self.gaps = np.array(self._g)

    # ---- no sorted/min/mean accessor exists on purpose ----
    def _gaps(self):
        d = np.diff(np.r_[self._x, self._x[0] + self.N])
        g = np.mod(d, self.N)
        g[g == 0] = self.N
        tot = g.sum()
        if abs(tot - self.N) > 1e-9:
            raise Lost(f"identity gaps sum to {tot}, not {self.N}: names crossed")
        return g

    def _match(self, p):
        """Bijective nearest-neighbour assignment, named wall -> detection."""
        d = np.abs(self._last[:, None] - p[None, :].astype(float))
        d = np.minimum(d, self.N - d)
        # cap must exceed the measured jitter and stay well under the minimum
        # wall separation, or a name could grab its neighbour. Rings here hold
        # walls 60+ blocks apart; the cap below is 13-18.
        cap = MAXV * self.every + JITTER
        out = np.full(self.n, -1)
        used = np.zeros(len(p), bool)
        order = np.dstack(np.unravel_index(np.argsort(d, axis=None), d.shape))[0]
        for i, j in order:
            if out[i] >= 0 or used[j] or d[i, j] > cap:
                continue
            out[i], used[j] = j, True
        if (out < 0).any():
            raise Lost(f"no detection within {cap:.0f} blocks for name(s) "
                       f"{np.flatnonzero(out < 0).tolist()}")
        return out

    def run(self, s, T):
        """Step to t0+T, or until a SHORTFALL in wall count persists. Self."""
        t = self._t0
        short_since = None
        for _ in range(T):
            s = step(s)
            t += 1
            if (t - self._t0) % self.every:
                continue
            ch, _ = walls(s, t, self.w)
            p = np.flatnonzero(ch)
            if len(p) < self.n:
                # candidate annihilation -- believe it only if it persists
                if short_since is None:
                    short_since = t
                if (t - short_since) >= self.persist * self.every:
                    self.t_event = short_since
                    break
                continue
            if short_since is not None:
                self.holes.append((short_since, t))
                self.blinks += 1
                short_since = None
            if len(p) > self.n:
                self.blinks += 1
            try:
                j = self._match(p)
            except Lost:
                # an unmatchable detection set is treated like a shortfall --
                # a candidate event that has to persist -- and counted apart
                # from blinks so a cap that is too tight shows up as a number
                # instead of as a short track.
                self.unmatched += 1
                if short_since is None:
                    short_since = t
                if (t - short_since) >= self.persist * self.every:
                    self.t_event = short_since
                    break
                continue
            new = p[j].astype(float)
            d = new - self._last
            d -= self.N * np.round(d / self.N)        # minimal circular step
            self._x = self._x + d
            self._last = new
            self._t.append(t)
            self._pos.append(self._x.copy())
            self._g.append(self._gaps())
        self.t = np.array(self._t)
        self.pos = np.array(self._pos)
        self.gaps = np.array(self._g)
        self.state = s
        return self


def pooled(tracks, frame):
    """Gaps of every track at a given frame index -- or a refusal.

    THE WHOLE POINT. If any track stopped before `frame`, this does not quietly
    average the ones that got there. It raises, with the counts, and the caller
    has to decide in the open. `nullwalk.py` measured what the quiet version
    does: a forceless walk reports (60, 90, 150) -> (80.5, 89.9, 129.6).
    """
    short = [k for k, tr in enumerate(tracks) if len(tr.t) <= frame]
    if short:
        raise Lost(f"{len(short)} of {len(tracks)} tracks never reached frame "
                   f"{frame}; averaging them would be a survivorship filter. "
                   f"Lengths: {sorted(len(tr.t) for tr in tracks)}")
    return np.array([tr.gaps[frame] for tr in tracks])


# ------------------------------------------------------------------ gates

def _ring3(N, cuts, offs, rng):
    s = locked_ring(N, rng)
    for i in range(len(cuts)):
        b = locked_ring(N, rng)
        for _ in range(offs[i]):
            b = step(b)
        lo, hi = cuts[i], cuts[(i + 1) % len(cuts)]
        if hi > lo:
            s[3 * lo:3 * hi] = b[3 * lo:3 * hi]
        else:
            s[3 * lo:], s[:3 * hi] = b[3 * lo:], b[:3 * hi]
    return s


def gate_sum(rng):
    """G1: identity gaps sum to N for 4000 steps on a live three-wall ring."""
    N = 300
    for offs in ((0, 2, 1), (0, 1, 2)):
        s = _ring3(N, [0, 60, 150], offs, rng)
        tr = Track(s, N, every=20).run(s, 4000)
        if len(tr.t) < 20:
            continue
        tot = tr.gaps.sum(1)
        if not np.allclose(tot, N):
            return f"gate_sum: gaps sum to {tot.min()}..{tot.max()}, not {N}"
    return None


def gate_velocity(rng):
    """G2: recover fire 341's isolated-wall velocities, which were measured by a
    DIFFERENT instrument (walls.py, phase-boundary crossing of a smoothed
    indicator) on a different construction (a two-domain splice). Charge 1 at
    -0.2414 and charge 2 at +0.1932, 74 sigma apart over 38 walls.

    The gate asks for the SIGN and for 15% on the magnitude. Tighter than that
    would be a gate about the smoothing width rather than about the physics.
    """
    N, T = 1200, 1500
    want = {1: -0.2414, 2: +0.1932}
    got = {}
    for off, ch_expect in ((1, 1), (2, 2)):
        vs = []
        for _ in range(6):
            s = locked_ring(N, rng)
            b = locked_ring(N, rng)
            for _ in range(off):
                b = step(b)
            s[3 * (N // 4):3 * (3 * N // 4)] = b[3 * (N // 4):3 * (3 * N // 4)]
            try:
                tr = Track(s, N, every=25).run(s, T)
            except Lost:
                continue
            if len(tr.t) < 20 or tr.n != 2:
                continue
            # the two walls of a splice carry charges q and 3-q; take the one
            # whose charge matches, by NAME, not by which moved further
            for i, q in enumerate(tr.charges):
                if q == ch_expect:
                    v = np.polyfit(tr.t[5:], tr.pos[5:, i], 1)[0]
                    vs.append(v)
        if not vs:
            return f"gate_velocity: no usable splice for charge {ch_expect}"
        got[ch_expect] = float(np.mean(vs))
    msg = "  ".join(f"q={q}: {v:+.4f} (want {want[q]:+.4f})" for q, v in got.items())
    for q, v in got.items():
        if np.sign(v) != np.sign(want[q]) or abs(abs(v) - abs(want[q])) > 0.15 * abs(want[q]):
            return f"gate_velocity: {msg}"
    return "OK " + msg


def gate_stop(rng):
    """G3: the hard stop fires, and for the RIGHT reason.

    THE FIRST VERSION OF THIS GATE PASSED FOR THE WRONG REASON and is the whole
    reason the gate is written this way. It put a mixed-charge pair 12 blocks
    apart and reported "stopped at t=5 with 1 clean frame". But `resolution.py`
    measured that this counter ERASES any domain narrower than w = 15 blocks, so
    a 12-block arc was never two resolved walls -- the stop fired on the
    counter's resolution limit, not on an annihilation. A gate that goes green
    for a reason other than the one it names is worse than no gate.

    So: build the arc at 60 blocks, four times the resolution; require that the
    track runs CLEANLY for a stretch first, and that the stop lands within a
    factor of two of the time the measured closing speed predicts. The pair
    closes at |v1| + |v2| = 0.434 blocks/step, so 60 blocks is ~138 steps.
    """
    N = 600
    pred = 60 / (0.2414 + 0.1932)
    for _ in range(10):
        s = _ring3(N, [0, 60], (0, 1), rng)
        try:
            tr = Track(s, N, every=5).run(s, 3000)
        except Lost:
            continue
        if tr.n != 2 or tr.t_event is None:
            continue
        if len(tr.t) < 8:
            continue
        if not (0.5 * pred < tr.t_event < 2.0 * pred):
            return (f"gate_stop: stopped at t={tr.t_event}, but the measured "
                    f"closing speed predicts {pred:.0f} +- a factor of two")
        return (f"OK  {len(tr.t)} clean frames, then stopped at t={tr.t_event} "
                f"(closing speed predicts {pred:.0f})")
    return "gate_stop: never saw a resolvable pair close inside 3000 steps"


def gate_blink(rng):
    """G4: a (1,1,1) ring survives the blinks. The pre-fix tracker died after a
    median of 650 steps on exactly this construction; flicker.py showed the same
    rings still carrying three walls at t = 6000."""
    N, T = 300, 6000
    lens = []
    for _ in range(6):
        s = _ring3(N, [0, 60, 150], (0, 1, 2), rng)
        try:
            tr = Track(s, N, every=5)
        except Lost:
            continue
        if tr.n != 3:
            continue
        tr.run(s, T)
        lens.append(tr.t[-1] - tr.t[0])
    if not lens:
        return "gate_blink: no (1,1,1) ring built"
    med = int(np.median(lens))
    if med < 3000:
        return (f"gate_blink: median tracked span only {med} steps "
                f"(spans {sorted(lens)}); the blinks are still stopping it")
    return (f"OK  median tracked span {med} of {T} steps over {len(lens)} rings "
            f"(pre-fix median was 650)")


if __name__ == "__main__":
    rng = np.random.default_rng(342342)
    why = gate_sum(rng)
    assert why is None, why
    print("G1 identity gaps sum to N for every frame of a live three-wall ring")
    for k, g in enumerate((gate_velocity, gate_stop, gate_blink), start=2):
        r = g(rng)
        assert r and r.startswith("OK"), r
        print(f"G{k} {r[3:]}")
    print("\ngates pass")
