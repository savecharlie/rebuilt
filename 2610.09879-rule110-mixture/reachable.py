#!/usr/bin/env python3
"""reachable.py -- how much of each rule table does the construction actually use?

arXiv:2610.09879 emulates rule 110 with (36, 37, 5) at residues (0, 1, 2) mod 3,
and its footnote 1 says rule 36 may be replaced by rule 54 "as the two rules only
differ in positions that are never visited".

That invites the exact question the paper does not ask: how many of the 256^3
rule triples emulate rule 110 under this encoding, and are the trivial ones
special among them?

METHOD. The mixture has radius 1 and the claim is about 3 steps, so every
(residue, neighbourhood) pair that can ever occur already occurs among the
encoded configurations of a few blocks. Enumerate them exhaustively, check the
visited set SATURATES in L, then ask of each visited pair whether its output is
FORCED (flip it; does the emulation break?). The count of working triples is then
2^(24 - #forced), and everything else follows by enumeration.
"""
import itertools

RES = (0, 1, 2)


def table(r):
    return {(k >> 2 & 1, k >> 1 & 1, k & 1): (r >> k) & 1 for k in range(8)}


def nbr(s, i):
    n = len(s)
    return (s[(i - 1) % n], s[i], s[(i + 1) % n])


def step(s, tabs):
    return tuple(tabs[i % 3][nbr(s, i)] for i in range(len(s)))


def enc(bits):
    out = []
    for b in bits:
        out += [0, 0, 0] if b else [1, 0, 0]
    return tuple(out)


def dec(blk):
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


def step110(s):
    t = table(110)
    return tuple(t[nbr(s, i)] for i in range(len(s)))


BASE = {0: table(36), 1: table(37), 2: table(5)}

print("1. which (residue, neighbourhood) pairs the construction ever visits")
prev = None
for L in range(2, 11):
    vis = set()
    for bits in itertools.product((0, 1), repeat=L):
        s = enc(bits)
        for _ in range(3):
            for i in range(len(s)):
                vis.add((i % 3, nbr(s, i)))
            s = step(s, BASE)
    print(f"   L={L:<3} visited pairs: {len(vis):>3}/24"
          + ("  (saturated)" if vis == prev else ""))
    if vis == prev:
        break
    prev = vis
VIS = prev
for r in RES:
    seen = sorted(k for rr, k in VIS if rr == r)
    print(f"   residue {r} (rule {(36,37,5)[r]}): uses {len(seen)}/8 entries: {seen}")

print("\n2. is every visited entry FORCED?  (flip one bit, does the claim break?)")
def emulates(tabs, Lmax=9):
    for L in range(1, Lmax + 1):
        for bits in itertools.product((0, 1), repeat=L):
            s = enc(bits)
            for _ in range(3):
                s = step(s, tabs)
            if dec(s) != step110(bits):
                return False
    return True


assert emulates(BASE), "baseline must emulate"
forced, free = [], []
for r in RES:
    for k in itertools.product((0, 1), repeat=3):
        tabs = {q: dict(BASE[q]) for q in RES}
        tabs[r][k] ^= 1
        (forced if not emulates(tabs) else free).append((r, k))
print(f"   forced entries: {len(forced)}   free entries: {len(free)}")
print(f"   forced == visited: {set(forced) == VIS}")
print(f"   free entries (never visited): "
      f"{sorted((r, ''.join(map(str, k))) for r, k in free)}")

print("\n3. so the number of rule triples that EXACTLY emulate rule 110 is 2^free")
nfree = len(free)
print(f"   2^{nfree} = {2 ** nfree}   out of 256^3 = {256**3}")
# enumerate them all and ask what they look like
triples = []
for assign in itertools.product((0, 1), repeat=nfree):
    tabs = {q: dict(BASE[q]) for q in RES}
    for (r, k), v in zip(free, assign):
        tabs[r][k] = v
    nums = tuple(sum((1 << (4 * k[0] + 2 * k[1] + k[2])) for k in tabs[r] if tabs[r][k])
                 for r in RES)
    triples.append(nums)
assert len(set(triples)) == 2 ** nfree
print(f"   every one verified to emulate: "
      f"{all(emulates({r: table(t[r]) for r in RES}, 7) for t in triples)}")
print(f"   the full list, by residue:")
for r in RES:
    print(f"     residue {r}: {sorted({t[r] for t in triples})}")


def symmetric(rule):
    t = table(rule)
    return all(t[(a, b, c)] == t[(c, b, a)] for a in (0, 1) for b in (0, 1) for c in (0, 1))


print("\n4. which of them are ALL mirror-symmetric (the paper's point)?")
sym = [t for t in triples if all(symmetric(x) for x in t)]
print(f"   {len(sym)} of {len(triples)}: {sym}")

print("\n5. triviality of the parts -- measured, per rule that appears")
def behaviour(rule, trials=400):
    import random
    random.seed(11)
    pers = set()
    for _ in range(trials):
        n = 3 * random.randint(3, 12)
        s = tuple(random.randint(0, 1) for _ in range(n))
        seen = {}
        for t_ in range(3000):
            if s in seen:
                pers.add(t_ - seen[s])
                break
            seen[s] = t_
            s = tuple(table(rule)[nbr(s, i)] for i in range(n))
    return max(pers)


allr = sorted({x for t in triples for x in t})
for rule in allr:
    print(f"   rule {rule:<4} mirror-symmetric {symmetric(rule)!s:<5} "
          f"longest cycle from random starts: {behaviour(rule)}")
