#!/usr/bin/env python3
"""sweep32b.py -- velocity AND core structure per implementation, done with a
tracker instead of a per-frame label.

WHY THIS REPLACES sweep32.py, WHICH IS RETRACTED. sweep32.py assigned each
frame's two detections to charge 1 and charge 2 by reading the charge off that
frame.  A blink (fire 342: ~2% of frames carry a spurious extra wall, and a
single wall can momentarily resolve as two of the other charge) therefore swaps
the two position series mid-run, and a least-squares slope through a swapped
series is a number about the swap.  It showed as error bars of +-0.07 on a
velocity whose seed-to-seed spread is +-0.006, and as charge-2 velocities of
+0.12, -0.09, -0.11 for walls that veltri.py measures at +0.18.  The symptom was
visible in the error bar; I would not have caught it from the mean.

So: track the two walls by IDENTITY with a jump cap of 20 blocks (jump.py
measured the detector's own frame-to-frame scatter at up to 8 blocks against a
wall speed of 0.25), take each tracked wall's charge and its two domain phases
as the MODE over the whole run, and drop a wall whose modal charge has less than
90% support.  Then the raw-label profile is accumulated against phases that do
not flicker.

*** THE `amb` CHANNEL CONFLATES TWO OPPOSITE THINGS. DO NOT READ IT AS ONE. ***
`labels()` returns -1 for the form `100`, which is in ALL THREE alphabets and so
has no exclusive beat, AND for a form in NO alphabet, which is illegal.  Under
{36,38} x {37} x {13,29} the wide region at the wall is 100% the ILLEGAL form
`001`; the baseline's few unlabelled blocks at a wall centre really are `100`.  I
read every -1 as the all-beats kind and wrote it into five files before
`invade.py` measured it.  Use invade.py for anything about that region; `amb` is
only a crude "something here is not a clean domain" channel.

GATE S -- four numbers measured in this directory by three earlier, independent
instruments must come back:
    v(charge 1) = -0.2414   v(charge 2) = +0.1932       (fire 341, charge_velocity.py)
    MID peak, charge 1 = 0.085 +- 0.007                 (fire 342, core.py)
    MID peak, charge 2 = 0.058 +- 0.003                 (fire 342, core.py)
"""
import json
import sys

import numpy as np

from tri import MAIN_CUBE, Kit

N, W, T, BURN, R, JUMP = 1800, 15, 800, 150, 45, 20
SEEDS = range(8)
REF = {"v1": -0.2414, "v2": +0.1932, "mid1": 0.085, "mid2": 0.058}


def one(kit, seed):
    rng = np.random.default_rng(4000 + 37 * seed)
    s = kit.two_domain(N, 1, rng)
    ch, _ = kit.walls(s, 0, W)
    p = np.flatnonzero(ch)
    if len(p) != 2:
        return None
    cur = p.astype(float)
    pos, chg, lf, rt, frames = [], [], [], [], []
    states = []
    for t in range(T):
        ch, _ = kit.walls(s, t, W)
        ph = kit.phase_field(s, t, W)
        q = np.flatnonzero(ch).astype(float)
        if len(q) == 0:
            return None
        new = cur.copy()
        for i in (0, 1):
            d = (q - cur[i] + N / 2) % N - N / 2
            j = int(np.abs(d).argmin())
            if abs(d[j]) <= JUMP:
                new[i] = cur[i] + d[j]
        cur = new
        j0 = (np.round(cur).astype(int)) % N
        pos.append(cur.copy())
        chg.append([int(ch[j0[i]]) for i in (0, 1)])
        lf.append([int(ph[(j0[i] - 40) % N]) for i in (0, 1)])
        rt.append([int(ph[(j0[i] + 40) % N]) for i in (0, 1)])
        frames.append(t)
        if t >= BURN:
            states.append((t, kit.labels(s, t), j0.copy()))
        s = kit.step(s)
    pos = np.array(pos); chg = np.array(chg); lf = np.array(lf); rt = np.array(rt)
    frames = np.array(frames)

    def mode(a):
        v, c = np.unique(a, return_counts=True)
        return int(v[c.argmax()]), c.max() / len(a)

    out = {}
    for i in (0, 1):
        c, sup = mode(chg[BURN:, i])
        l, supl = mode(lf[BURN:, i])
        r, supr = mode(rt[BURN:, i])
        if sup < 0.9 or supl < 0.9 or supr < 0.9 or c not in (1, 2):
            continue
        m = (l + 2 * c) % 3              # the third phase, neither domain's
        if (r - l) % 3 != c:
            continue                     # the modal phases must agree with the modal charge
        v = np.polyfit(frames[BURN:], pos[BURN:, i], 1)[0]
        mid = np.zeros(2 * R + 1); amb = np.zeros(2 * R + 1); n = 0
        for (t, lab, j0) in states:
            idx = (np.arange(-R, R + 1) + j0[i]) % N
            L = lab[idx]
            mid += (L == m); amb += (L < 0); n += 1
        out[c] = (float(v), mid / n, amb / n)
    return out or None


def run(kit, seeds=SEEDS):
    acc = {1: [], 2: []}
    drop = 0
    for sd in seeds:
        r = one(kit, sd)
        if r is None:
            drop += 1
            continue
        for c, val in r.items():
            acc[c].append(val)
    e = {"dropped": drop}
    for c in (1, 2):
        if not acc[c]:
            e[f"v{c}"] = None
            continue
        v = np.array([a[0] for a in acc[c]])
        M = np.array([a[1] for a in acc[c]])
        A = np.array([a[2] for a in acc[c]])
        e[f"v{c}"] = [float(v.mean()), float(v.std(ddof=1) / np.sqrt(len(v)))
                      if len(v) > 1 else 0.0]
        e[f"mid{c}"] = M.mean(0).tolist()
        e[f"mid{c}_se"] = (M.std(0, ddof=1) / np.sqrt(len(M))).tolist() if len(M) > 1 \
            else [0.0] * M.shape[1]
        e[f"amb{c}"] = A.mean(0).tolist()
        e[f"n{c}"] = int(len(v))
    return e


if __name__ == "__main__":
    print("GATE S -- baseline (36,37,5) against four numbers from three earlier instruments")
    b = run(Kit((36, 37, 5)))
    ok = True
    for c in (1, 2):
        v, ve = b[f"v{c}"]
        pk = max(b[f"mid{c}"]); pki = int(np.argmax(b[f"mid{c}"])) - R
        gv = abs(v - REF[f"v{c}"]) < 0.02
        gm = abs(pk - REF[f"mid{c}"]) < 0.03
        ok &= gv and gm
        print(f"   charge {c}: v {v:+.4f}+-{ve:.4f} (ref {REF[f'v{c}']:+.4f}) {'ok' if gv else 'FAIL'}"
              f"   MID peak {pk:.4f} at offset {pki:+d} (ref {REF[f'mid{c}']:.3f}) {'ok' if gm else 'FAIL'}"
              f"   [{b[f'n{c}']} walls, {b['dropped']} seeds dropped]")
    print(f"   -> {'PASS' if ok else 'FAIL'}")
    if not ok:
        sys.exit(1)

    res = {}
    for i, tri in enumerate(MAIN_CUBE):
        e = run(Kit(tri))
        res[str(tri)] = e
        def fm(c):
            if e.get(f"v{c}") is None:
                return "   none  "
            return f"{e[f'v{c}'][0]:+.4f}+-{e[f'v{c}'][1]:.4f} pk{max(e[f'mid{c}']):.3f}"
        print(f"[{i+1}/32] {str(tri):<16} c1 {fm(1)}   c2 {fm(2)}   drop={e['dropped']}",
              flush=True)
    json.dump({"R": R, "N": N, "T": T, "BURN": BURN, "data": res},
              open("sweep32b.json", "w"))
    print("wrote sweep32b.json")
