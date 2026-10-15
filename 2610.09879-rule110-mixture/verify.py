#!/usr/bin/env python3
"""verify.py -- rebuild arXiv:2610.09879 (Noronha, Oct 7 2026) from its own words.

CLAIM. Assign ECA rule 36 to sites i%3==0, rule 37 to i%3==1, rule 5 to i%3==2.
All three are mirror-symmetric and individually trivial. Encode logical 0 as the
block [1,0,0] and logical 1 as [0,0,0]. Then THREE steps of the mixture equal ONE
step of rule 110 -- so the mixture is Turing complete, and the chirality of 110
comes from the ARRANGEMENT, not from any rule.

This is a claim with no free parameters and no tolerance: it is either exactly
true on every configuration or it is false. So the check is exhaustive, in
integer arithmetic, with the instrument gated against a known answer first.

GATES (each must pass before the next claim is tested):
  G1  my ECA stepper reproduces rule 110's published rule table and its
      left-moving glider background -- the stepper is validated before use.
  G2  rules 36, 37, 5 are mirror-symmetric; 110 is not.
  G3  rules 36, 37, 5 are individually trivial (every orbit is eventually
      periodic with a SHORT period, from random starts).
Then:
  C1  the emulation, exhaustively over every configuration of L = 1..14 logical
      cells (periodic), plus random large L.
  C2  the footnote: rule 36 may be replaced by rule 54 because they differ only
      at neighbourhoods never visited at i%3==0 sites.
  C3  the footnote's combinatorics: no period-<=5 two-symbol arrangement is
      chiral, (a,a,b,a,b,b) at period 6 is, and (a,b,c) at period 3 is.
"""
import itertools
import random

random.seed(20261013)


def table(rule):
    """rule number -> dict (l,c,r) -> out, Wolfram convention k = 4l+2c+r."""
    return {(k >> 2 & 1, k >> 1 & 1, k & 1): (rule >> k) & 1 for k in range(8)}


T = {r: table(r) for r in (5, 36, 37, 54, 110, 118)}


def step_uniform(state, rule):
    n = len(state)
    t = T[rule]
    return tuple(t[(state[(i - 1) % n], state[i], state[(i + 1) % n])] for i in range(n))


def step_mixture(state, rules=(36, 37, 5)):
    """rules[i % 3] applied at site i."""
    n = len(state)
    assert n % 3 == 0, "the ring length must be a multiple of 3 for the tiling to close"
    return tuple(T[rules[i % 3]][(state[(i - 1) % n], state[i], state[(i + 1) % n])]
                 for i in range(n))


def encode(bits):
    """logical 0 -> [1,0,0]; logical 1 -> [0,0,0]."""
    out = []
    for b in bits:
        out += [0, 0, 0] if b else [1, 0, 0]
    return tuple(out)


def decode(block):
    """inverse of encode; raises if a block is neither [100] nor [000]."""
    out = []
    for j in range(0, len(block), 3):
        b = block[j:j + 3]
        if b == (1, 0, 0):
            out.append(0)
        elif b == (0, 0, 0):
            out.append(1)
        else:
            raise ValueError(f"block {j//3} is {b}, not a legal codeword")
    return tuple(out)


fails = []

# ---------------------------------------------------------------- G1
print("G1  the stepper, against THREE published invariants of rule 110")
pub110 = {(1,1,1):0, (1,1,0):1, (1,0,1):1, (1,0,0):0,
          (0,1,1):1, (0,1,0):1, (0,0,1):1, (0,0,0):0}
g1a = all(T[110][k] == v for k, v in pub110.items())
print(f"    (a) rule 110 table == 01101110 = 110: {g1a}")


def mirror_rule(rule):
    t = table(rule)
    return sum((1 << k) for k in range(8) if t[(k & 1, k >> 1 & 1, k >> 2 & 1)])


# the documented mirror-conjugate pair 110 <-> 124, and the constituents being
# their OWN mirrors -- the symmetry claim reached a second, independent way
g1b = mirror_rule(110) == 124 and mirror_rule(124) == 110
g1c = [mirror_rule(r) for r in (36, 37, 5, 54)] == [36, 37, 5, 54]
print(f"    (b) mirror(110) == 124 and mirror(124) == 110: {g1b}")
print(f"        mirror(36,37,5,54) == themselves: {g1c}")
# Cook's ether: rule 110's background has spatial period 14, temporal period 7.
# Search EVERY 14-ring for step^7 == some cyclic shift.  Nothing is told to it.
ether = [b for b in itertools.product((0, 1), repeat=14)
         if 0 < sum(b) < 14
         and any(_s == tuple(b[(i + k) % 14] for i in range(14)) for k in range(14)
                 for _s in [__import__("functools").reduce(
                     lambda a, _: step_uniform(a, 110), range(7), b)])]
g1d = (len(ether) == 14 and len({sum(b) for b in ether}) == 1
       and sum(ether[0]) == 8)
print(f"    (c) nontrivial period-7 patterns on a 14-ring: {len(ether)} "
      f"= one orbit of 14 shifts at density {sum(ether[0])}/14: {g1d}")
print(f"        the ether: {''.join('#' if v else '.' for v in ether[0])}")
if not (g1a and g1b and g1c and g1d):
    fails.append("G1")

# ---------------------------------------------------------------- G2
print("\nG2  mirror symmetry")
def symmetric(rule):
    t = T[rule]
    return all(t[(l, c, r)] == t[(r, c, l)] for l in (0, 1) for c in (0, 1) for r in (0, 1))
for r in (5, 36, 37, 54):
    ok = symmetric(r)
    print(f"    rule {r:<3} mirror-symmetric: {ok}")
    if not ok:
        fails.append(f"G2/{r}")
for r in (110, 118):
    ok = symmetric(r)
    print(f"    rule {r:<3} mirror-symmetric: {ok}   (must be False)")
    if ok:
        fails.append(f"G2/{r}")
# and the structural reading: all four are (l==r) AND g(l,c)
asym = [(r, l, c, rr) for r in (5, 36, 37, 54)
        for l in (0, 1) for c in (0, 1) for rr in (0, 1)
        if l != rr and T[r][(l, c, rr)] != 0]
print(f"    all four output 0 whenever left != right, so each is "
      f"(l==r) AND g(l,c): {not asym}")

# ---------------------------------------------------------------- G3
print("\nG3  each constituent rule is individually trivial (short cycles from random starts)")
for r in (36, 37, 5):
    worst = 0
    for _ in range(300):
        n = 3 * random.randint(4, 20)
        s = tuple(random.randint(0, 1) for _ in range(n))
        seen = {}
        for t_ in range(2000):
            if s in seen:
                worst = max(worst, t_ - seen[s])
                break
            seen[s] = t_
            s = step_uniform(s, r)
    print(f"    rule {r:<3} longest cycle over 300 random rings: {worst}")
    if worst > 4:
        fails.append(f"G3/{r}")

# ---------------------------------------------------------------- C1
print("\nC1  THE CLAIM: 3 mixture steps == 1 rule-110 step, exhaustively")
bad = []
for L in range(1, 15):
    for bits in itertools.product((0, 1), repeat=L):
        want = step_uniform(bits, 110)
        blk = encode(bits)
        for _ in range(3):
            blk = step_mixture(blk)
        try:
            got = decode(blk)
        except ValueError as exc:
            bad.append((L, bits, str(exc)))
            continue
        if got != want:
            bad.append((L, bits, f"{got} != {want}"))
    tested = 2 ** L
    print(f"    L={L:<3} {tested:>6} configurations   mismatches so far: {len(bad)}")
    if bad:
        break
if not bad:
    for L in (50, 101, 257, 1000):
        bits = tuple(random.randint(0, 1) for _ in range(L))
        blk = encode(bits)
        for _ in range(3):
            blk = step_mixture(blk)
        ok = decode(blk) == step_uniform(bits, 110)
        print(f"    random L={L:<5} {ok}")
        if not ok:
            bad.append((L, "random", ""))
    # and many steps in a row, to be sure it is not a one-step coincidence
    bits = tuple(random.randint(0, 1) for _ in range(300))
    blk, ref = encode(bits), bits
    for t_ in range(200):
        for _ in range(3):
            blk = step_mixture(blk)
        ref = step_uniform(ref, 110)
        if decode(blk) != ref:
            bad.append(("iterate", t_, ""))
            break
    print(f"    200 consecutive rule-110 steps on L=300 tracked exactly: {not bad}")
print(f"    C1 {'CONFIRMED' if not bad else f'FAILED {bad[:3]}'}")
if bad:
    fails.append("C1")

# ---------------------------------------------------------------- C2
print("\nC2  footnote: rule 36 -> rule 54 (they differ only where never visited)")
diff = [k for k in T[36] if T[36][k] != T[54][k]]
print(f"    36 and 54 differ at neighbourhoods: {sorted(diff)}")
visited = set()
for _ in range(400):
    L = random.randint(4, 40)
    bits = tuple(random.randint(0, 1) for _ in range(L))
    blk = encode(bits)
    for _ in range(3 * 40):
        n = len(blk)
        for i in range(0, n, 3):                      # only the i%3==0 sites
            visited.add((blk[(i - 1) % n], blk[i], blk[(i + 1) % n]))
        blk = step_mixture(blk)
print(f"    neighbourhoods ever seen at an i%3==0 site: {sorted(visited)}")
c2 = not (set(diff) & visited)
print(f"    the differing neighbourhoods are never visited: {c2}")
# and the strong form: swap the rule and re-run the whole claim
bad54 = []
for L in range(1, 13):
    for bits in itertools.product((0, 1), repeat=L):
        blk = encode(bits)
        for _ in range(3):
            blk = step_mixture(blk, rules=(54, 37, 5))
        if decode(blk) != step_uniform(bits, 110):
            bad54.append((L, bits))
print(f"    with rule 54 in place of 36, the emulation still holds for L<=12: "
      f"{not bad54}")
if not c2 or bad54:
    fails.append("C2")

# ---------------------------------------------------------------- C3
print("\nC3  footnote combinatorics: the shortest CHIRAL periodic arrangement")
def chiral(seq):
    """no cyclic shift of the reversal equals the sequence."""
    n = len(seq)
    rev = seq[::-1]
    return all(tuple(rev[(i + s) % n] for i in range(n)) != tuple(seq) for s in range(n))
for k in (2, 3):
    first = None
    for n in range(1, 9):
        hits = [s for s in itertools.product(range(k), repeat=n) if chiral(s)]
        if hits and first is None:
            first = (n, hits[0], len(hits))
        print(f"    {k} symbols, period {n}: {len(hits):>5} chiral arrangements")
        if first and n >= first[0] + 1:
            break
    print(f"    -> {k} symbols: shortest chiral period = {first[0]}, e.g. {first[1]}")
c3 = chiral((0, 0, 1, 0, 1, 1)) and not any(
    chiral(s) for n in range(1, 6) for s in itertools.product(range(2), repeat=n))
print(f"    (a,a,b,a,b,b) is chiral AND no 2-symbol period <= 5 is: {c3}")
print(f"    (a,b,c) at period 3 is chiral: {chiral((0,1,2))}")
if not c3:
    fails.append("C3")

print("\n" + "=" * 62)
print("ALL GATES AND CLAIMS PASS" if not fails else f"FAILURES: {fails}")
