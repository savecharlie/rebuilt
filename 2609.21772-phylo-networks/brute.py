"""Count rooted binary phylogenetic networks by brute force — no component graphs.

Independent instrument for Yu & Zhang, arXiv:2609.21772. Their count comes from
classifying 79 component graphs into ten groups and summing lemmas. This one
knows nothing about components: it grows every network top-down from the root,
one node at a time, and throws away duplicates with nauty.

A partial network has *stubs*: edges whose head is not decided yet. A step takes
one stub and makes it
  - a leaf with an unused label,
  - a tree node (two new stubs below it),
  - a reticulation, by merging it with another stub whose tail is different
    (same tail would be a parallel edge, which the definition forbids),
    with one new stub below.
Every network is reached (build it in any topological order), and nauty's
canonical certificate, with the leaf labels as singleton colour classes,
collapses the orderings. Acyclicity is free: new nodes only ever point down.

Validate against a number already known before believing a new one
(the rebuilding skill). Known independently of the authors:
  k=0 : (2n-3)!!           trees
  k=1 : n(2n)!/(2^n n!) - 2^(n-1) n!      (Zhang 2019)
  k=2 : Mansouri et al 2020 formula

Iris, Sep 21 2026.
"""
import sys
from math import factorial
from pynauty import Graph, certificate

R, T, H, S = 'R', 'T', 'H', 'S'


def canon(types, out):
    """types: list of type tags (int for leaf label). out: list of child lists."""
    n = len(types)
    classes = {}
    for v, t in enumerate(types):
        classes.setdefault(t, set()).add(v)
    order = [R, T, H, S] + sorted(k for k in classes if isinstance(k, int))
    col = [classes[k] for k in order if k in classes]
    sig = tuple((k, len(classes[k])) for k in order if k in classes)
    g = Graph(n, directed=True, adjacency_dict={v: list(out[v]) for v in range(n)},
              vertex_coloring=col)
    return sig, certificate(g)


def count(n, k):
    # state: (types tuple, out tuple of tuples, parent-of-stub map implicit)
    types = [R, S]
    out = [[1], []]
    start = (tuple(types), tuple(tuple(c) for c in out))
    level = {canon(types, out): start}
    need_T = n + k - 1
    finals = set()
    while level:
        nxt = {}
        for types, out in level.values():
            types = list(types)
            out = [list(c) for c in out]
            nT = types.count(T)
            nH = types.count(H)
            used = {t for t in types if isinstance(t, int)}
            stubs = [v for v, t in enumerate(types) if t == S]
            if not stubs:
                if nT == need_T and nH == k and len(used) == n:
                    finals.add(canon(types, out))
                continue
            parent = {}
            for u, cs in enumerate(out):
                for c in cs:
                    if types[c] == S:
                        parent[c] = u
            cand = []
            for s in stubs:
                # leaf
                for lab in range(1, n + 1):
                    if lab not in used:
                        t2 = types[:]; t2[s] = lab
                        cand.append((t2, [c[:] for c in out]))
                # tree node
                if nT < need_T:
                    t2 = types[:] + [S, S]; t2[s] = T
                    o2 = [c[:] for c in out] + [[], []]
                    o2[s] = [len(types), len(types) + 1]
                    cand.append((t2, o2))
                # reticulation: merge s with s2 (s2 > s to halve work)
                if nH < k:
                    for s2 in stubs:
                        if s2 <= s or parent[s2] == parent[s]:
                            continue
                        o2 = [c[:] for c in out]
                        p2 = parent[s2]
                        o2[p2] = [s if c == s2 else c for c in o2[p2]]
                        t2 = types[:]; t2[s] = H
                        # delete node s2: relabel everything above it down by one
                        del t2[s2]; del o2[s2]
                        o2 = [[c - (c > s2) for c in cs] for cs in o2]
                        sn = s - (s > s2)
                        t2.append(S); o2.append([])
                        o2[sn] = [len(t2) - 1]
                        cand.append((t2, o2))
            for t2, o2 in cand:
                remT = need_T - t2.count(T)
                remH = k - t2.count(H)
                if remT < 0 or remH < 0:
                    continue
                key = canon(t2, o2)
                if key not in nxt:
                    nxt[key] = (tuple(t2), tuple(tuple(c) for c in o2))
        level = nxt
    return len(finals)


def dfact(m):
    r = 1
    while m > 1:
        r *= m; m -= 2
    return r


def known(n, k):
    if k == 0:
        return dfact(2 * n - 3) if n > 1 else 1
    if k == 1:
        return n * factorial(2 * n) // (2 ** n * factorial(n)) - 2 ** (n - 1) * factorial(n)
    if k == 2:
        a = factorial(2 * n - 2) * (6 * n**4 + 19 * n**3 + 18 * n**2 - 4 * n - 6)
        b = 3 * 2 ** (n - 1) * factorial(n - 1)
        assert a % b == 0
        return a // b - 2 ** (n - 1) * factorial(n + 1) * (2 * n + 3)
    return None


if __name__ == '__main__':
    n, k = int(sys.argv[1]), int(sys.argv[2])
    c = count(n, k)
    print(f"n={n} k={k}  brute={c}  known={known(n, k)}", flush=True)
