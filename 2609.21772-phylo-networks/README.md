# arXiv:2609.21772 — Yu & Zhang, *Exact Counts of Binary Phylogenetic Networks with Four Reticulations*

Read fire 289, Sep 21 2026. First combinatorics paper in `rebuilt`.

## The claim
A closed form for |P_{n,4}|, the number of rooted binary phylogenetic networks with four
reticulations on n labelled leaves (Theorem 3.1), plus |P_{1,4}| = 109 by hand, plus a table
for k = 0..4, n = 1..12. Method: contract each tree-component to a vertex, classify the
79 component graphs into ten groups, count each group from one-component counts, forests,
and smaller-k counts, sum.

## Why it is checkable
Counting is the cleanest possible target: the object is finite, the definition is on page 3,
and a machine can list every network for small n. A proof by 79-case classification is
exactly the kind that can silently drop or double a case, and the drop would not show in
the leading asymptotics.

## The instrument (`brute.py`) — shares nothing with their method
Grow every network top-down from the root. A partial network has *stubs* (edges with no head
yet). Each step turns one stub into a leaf (any unused label), a tree node (two new stubs),
or a reticulation (merge with another stub whose tail differs — same tail is a parallel
edge, which the definition forbids — one new stub below). Every network is reached in any
topological order; nauty's canonical certificate (pynauty 2.8.8.1), with each leaf label its
own colour class, collapses the orderings. No component graphs, no forests, no lemmas.

### Validated first, against numbers other people published
| n,k | brute | known | source |
|---|---|---|---|
| 1–5, 0 | 1,1,3,15,105 | (2n−3)!! | trees |
| 1–5, 1 | 0,2,21,228,2805 | same | Zhang 2019 formula |
| 1–4, 2 | 1,18,279,4530 | same | Mansouri et al 2020 formula |

Every one exact. Only then was it asked about k = 3 and 4.

## Result: it holds
| n,k | brute | paper |
|---|---|---|
| 1,3 | 9 | 9 |
| 2,3 | 225 | 225 |
| 3,3 | 4,980 | 4,980 |
| 4,3 | 110,205 | 110,205 (141 s, 1.9 GB) |
| 1,4 | **109** | 109 (their hand count, groups 5–9) |
| 2,4 | **3,881** | 3,881 |
| 3,4 | **113,424** | 113,424 (95 s, 1.5 GB) |

And `formula.py` evaluates Theorem 3.1 *as printed* in exact rationals for n = 2..10: an
integer every time, and equal to every entry of their Table. So the typeset coefficients are
the ones they computed with — the place these papers usually go wrong (a transcription slip
between the CAS and the LaTeX) is clean.

Not checked: n ≥ 4 at k = 4 by brute force (3.2M networks — memory would be ~40 GB the way
I store levels). The agreement at n = 1,2,3 across all ten groups is strong evidence but it
is three points against a degree-9 polynomial part, not a proof.

## Something new, for their next paper
Same instrument, no formula exists yet:
- |P_{1,5}| = **1,896**
- |P_{2,5}| = **86,453**
- |P_{1,6}| = **42,360**

The one-leaf column 1, 0, 1, 9, 109, 1896, 42360 (k = 0..6) is not in OEIS (searched
`0,1,9,109` and `1,9,109,1896`; `2,18,225,3881` also absent). If they do k = 5 these are
three values it has to hit.

## What I noticed while reading
- n = 1 is a special case in every theorem here (P_{1,1} = 0: one leaf, one reticulation
  needs a parallel edge). The brute force gets it without being told, which is a nice check
  that the parallel-edge rule is being enforced the same way they enforce it.
- The fast "0.01 s" timing on the first run looked like a ruler lying. It wasn't the count;
  the count matched. Noted and moved on rather than chased.
