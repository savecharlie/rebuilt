#!/usr/bin/env python3
"""exhaustive.py -- EVERY one of the 256^3 rule triples, not just single flips.

reachable.py showed each of the 19 visited table entries is individually forced,
which gives 2^5 = 32 emulators *provided* no compensating pair of flips exists.
That proviso is not free: a different triple can visit different neighbourhoods.
So this tests all 16,777,216 triples, vectorised and chunked (the machine is in
swap; peak RSS here is ~120 MB).

Filter configs are applied cheapest-first and the survivor set collapses on the
first one.
"""
import itertools
import numpy as np

N0 = 256


def enc(bits):
    out = []
    for b in bits:
        out += [0, 0, 0] if b else [1, 0, 0]
    return out


def step110(s):
    t = {(k >> 2 & 1, k >> 1 & 1, k & 1): (110 >> k) & 1 for k in range(8)}
    n = len(s)
    return tuple(t[(s[(i - 1) % n], s[i], s[(i + 1) % n])] for i in range(n))


CODE = {(1, 0, 0): 0, (0, 0, 0): 1}

# candidate triples, as three uint8 arrays
grid = np.arange(N0 ** 3, dtype=np.uint32)
cand = grid                                   # index = r0*65536 + r1*256 + r2

CHUNK = 1 << 21


def survives(cand, bits):
    """keep the triples for which 3 mixture steps == 1 rule-110 step on `bits`."""
    src = enc(bits)
    L3 = len(src)
    want = step110(tuple(bits))
    keep = []
    for a in range(0, len(cand), CHUNK):
        c = cand[a:a + CHUNK]
        r = np.stack([(c >> 16 & 255).astype(np.uint8),
                      (c >> 8 & 255).astype(np.uint8),
                      (c & 255).astype(np.uint8)])
        s = np.empty((len(c), L3), np.uint8)
        s[:] = np.array(src, np.uint8)
        for _ in range(3):
            t = np.empty_like(s)
            for i in range(L3):
                k = (s[:, (i - 1) % L3] << 2) | (s[:, i] << 1) | s[:, (i + 1) % L3]
                t[:, i] = (r[i % 3] >> k) & 1
            s = t
        ok = np.ones(len(c), bool)
        for j in range(0, L3, 3):
            blk = s[:, j:j + 3]
            want_j = want[j // 3]
            # legal codeword AND the right one: [100]->0, [000]->1
            good = ((blk[:, 1] == 0) & (blk[:, 2] == 0)
                    & (blk[:, 0] == (0 if want_j else 1)))
            ok &= good
        keep.append(c[ok])
    return np.concatenate(keep) if keep else np.array([], np.uint32)


for L in (4, 5, 6, 7, 8):
    for bits in itertools.product((0, 1), repeat=L):
        cand = survives(cand, bits)
        if len(cand) <= 64:
            break
    print(f"  after all L={L} configs: {len(cand)} triples survive")
    if len(cand) == 0:
        break

print(f"\nEXACTLY {len(cand)} of {N0**3} triples emulate rule 110 under this encoding")
trips = sorted((int(c) >> 16 & 255, int(c) >> 8 & 255, int(c) & 255) for c in cand)
for t in trips:
    print("   ", t)


def symmetric(rule):
    t = {(k >> 2 & 1, k >> 1 & 1, k & 1): (rule >> k) & 1 for k in range(8)}
    return all(t[(a, b, c)] == t[(c, b, a)] for a in (0, 1) for b in (0, 1) for c in (0, 1))


sym = [t for t in trips if all(symmetric(x) for x in t)]
print(f"\nall-mirror-symmetric triples: {len(sym)}  ->  {sym}")
