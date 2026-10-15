#!/usr/bin/env python3
"""invade.py -- the band is an ILLEGAL phase, and it eats the computation.

THE ERROR THIS FILE EXISTS TO CORRECT, written first because it reached five
other files before I caught it.  `labels()` returns -1 for TWO different reasons:
the block form `100`, which is in ALL THREE alphabets and so has no exclusive
beat, and a form in NO alphabet, which is simply illegal.  I read every -1 as the
first and wrote "a band of the form `100`, legal on all three beats and therefore
belonging to none" into NOTES.md, SECTION5_RESULTS.md, a letter, an essay and a
picture caption.

Measured instead: the wide band in (36,37,13) and (38,37,29) is **100% the form
`001`**, which is in no alphabet at all.  The baseline's small 5-block runs ARE
pure `100`.  The two are opposite in meaning and identical through the label
field, which is what a sentinel with two causes does.

So these implementations do not merely pin their walls.  **The defect manufactures
an illegal phase that is not a codeword at any beat, and it invades.**  A ring
containing one `001` block cannot be locked by definition, which is a complete
mechanism for lock32.py's zeros and is nothing to do with pinning.

MEASURED HERE: the growth law of the illegal region, its invasion speed, and
whether the all-illegal state is absorbing.
"""
import numpy as np

from tri import Kit, BASE_ALPH, MAIN_CUBE

ILLEGAL = [c for c in range(8)
           if not any(((c >> 2) & 1, (c >> 1) & 1, c & 1) in A
                      for A in BASE_ALPH.values())]
ALLBEATS = [4 * x + 2 * y + z for (x, y, z) in
            set.intersection(*[set(v) for v in BASE_ALPH.values()])]
N, T = 1200, 3000


def codes(s):
    b = s.reshape(-1, 3).astype(np.int64)
    return 4 * b[:, 0] + 2 * b[:, 1] + b[:, 2]


def frac_illegal(kit, seed, N=N, T=T, every=25):
    s = kit.two_domain(N, 1, np.random.default_rng(seed))
    out = []
    for t in range(T):
        if t % every == 0:
            c = codes(s)
            out.append((t, float(np.isin(c, ILLEGAL).mean()),
                        float(np.isin(c, ALLBEATS).mean())))
        s = kit.step(s)
    return np.array(out)


if __name__ == "__main__":
    print(f"form in NO alphabet (illegal): {['%d%d%d' % ((c>>2)&1,(c>>1)&1,c&1) for c in ILLEGAL]}")
    print(f"form in ALL alphabets        : {['%d%d%d' % ((c>>2)&1,(c>>1)&1,c&1) for c in ALLBEATS]}")

    print("\nGATE -- a locked ring with NO wall must never produce an illegal block,")
    print("under any of the 32 (the alphabets are closed under the dynamics).")
    worst = 0.0
    for tri in MAIN_CUBE:
        k = Kit(tri)
        s = k.locked_ring(600, np.random.default_rng(99))
        for _ in range(400):
            s = k.step(s)
            worst = max(worst, float(np.isin(codes(s), ILLEGAL).mean()))
    print(f"   max illegal fraction over 32 triples x 400 steps: {worst:.6f}"
          f"   -> {'PASS' if worst == 0 else 'FAIL'}")

    print(f"\nillegal-phase fraction of the ring over time (N={N} blocks, 3 seeds)")
    for tri in [(36, 37, 5), (36, 37, 13), (36, 37, 29), (38, 37, 13),
                (38, 37, 29), (52, 37, 29), (54, 37, 13)]:
        k = Kit(tri)
        runs = [frac_illegal(k, 31 + sd) for sd in range(3)]
        A = np.stack([r[:, 1] for r in runs])
        t = runs[0][:, 0]
        pick = [np.argmin(abs(t - x)) for x in (100, 500, 1000, 2000, 2975)]
        # invasion speed from the linear part, in BLOCKS per step
        m = (t >= 500) & (t <= 2000)
        sl = np.polyfit(t[m], A.mean(0)[m] * N, 1)[0]
        print(f"   {str(tri):<16} illegal frac at t=100/500/1000/2000/2975: "
              + " / ".join(f"{A.mean(0)[i]:.3f}" for i in pick)
              + f"   growth {sl:+.4f} blocks/step")

    print("\nis the all-illegal ring absorbing?  (every block `001`, 200 steps)")
    for tri in [(36, 37, 13), (38, 37, 29), (36, 37, 5)]:
        k = Kit(tri)
        s = np.tile(np.array([0, 0, 1], np.uint8), 400)
        same = True
        for _ in range(200):
            s2 = k.step(s)
            if not np.array_equal(s2, s):
                same = False
                break
            s = s2
        c = codes(s)
        print(f"   {str(tri):<16} fixed point: {same}   "
              f"illegal fraction after: {np.isin(c, ILLEGAL).mean():.3f}")


def classify_all(N=900, T=1500, seeds=3):
    """Every one of the 32: does its defect manufacture the illegal phase, and
    how fast does it spread?"""
    rows = []
    for tri in MAIN_CUBE:
        k = Kit(tri)
        A = []
        for sd in range(seeds):
            s = k.two_domain(N, 1, np.random.default_rng(31 + sd))
            f = []
            for t in range(T):
                if t % 25 == 0:
                    f.append(float(np.isin(codes(s), ILLEGAL).mean()))
                s = k.step(s)
            A.append(f)
        A = np.array(A)
        t = np.arange(A.shape[1]) * 25
        m = (t >= 300) & (t <= 1200)
        sl = np.polyfit(t[m], A.mean(0)[m] * N, 1)[0]
        rows.append((tri, float(A.mean(0)[-1]), float(sl),
                     float(A.mean(0)[-1] - A.mean(0)[len(t)//2])))
    return rows
