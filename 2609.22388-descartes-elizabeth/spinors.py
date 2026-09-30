"""The half of arXiv:2609.22388 that it asserts and never shows: what a spinor
is FOR in Descartes' theorem.  Exact integer / Gaussian-integer arithmetic.
LAW 1: every claim is gated against an independently known answer first.
"""
import sympy as sp
from math import isqrt

# ================================================================= CLAIM 1
J = sp.ones(4, 4)
Q = sp.eye(4) - sp.Rational(1, 2)*J
ev = Q.eigenvals()
sig = (sum(v for k, v in ev.items() if k > 0), sum(v for k, v in ev.items() if k < 0))
print("Q2 = I - (1/2)J  eigenvalues", {sp.nsimplify(k): v for k, v in ev.items()},
      "-> signature", sig, "= R^(3,1), Minkowski")
assert sig == (3, 1)
bs = sp.Matrix(sp.symbols('b1 b2 b3 b4'))
assert sp.expand((bs.T*Q*bs)[0] - (sum(x**2 for x in bs) - sp.Rational(1,2)*sum(bs)**2)) == 0
print("   => Descartes' theorem IS 'the curvature 4-vector is null'.")
print()

def gmul(a, b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def gadd(a, b): return (a[0]+b[0], a[1]+b[1])
def gsub(a, b): return (a[0]-b[0], a[1]-b[1])
def gsmul(k, a): return (k*a[0], k*a[1])
def gnorm(a): return a[0]*a[0] + a[1]*a[1]

def descartes_ok(q):
    b = [d[0] for d in q]; w = [d[1] for d in q]
    if sum(b)**2 != 2*sum(x*x for x in b): return False
    S = (0,0); T = (0,0)
    for x in w:
        S = gadd(S, x); T = gadd(T, gmul(x, x))
    return gmul(S, S) == gsmul(2, T)

def swap(q, i):
    sb = sum(d[0] for d in q) - q[i][0]
    sw = (0,0)
    for j, d in enumerate(q):
        if j != i: sw = gadd(sw, d[1])
    return (2*sb - q[i][0], gsub(gsmul(2, sw), q[i][1]))

seed = ((-1,(0,0)), (2,(-1,0)), (2,(1,0)), (3,(0,2)))
assert descartes_ok(seed), "GATE: seed is not a Descartes configuration"
print("GATE  seed (-1,2,2,3) satisfies real AND complex Descartes   OK")

seen = {tuple(sorted(seed))}
disks, quads, frontier = {}, [], [seed]
for _ in range(6):
    nxt = []
    for q in frontier:
        assert descartes_ok(q), "a generated quadruple broke Descartes"
        quads.append(q)
        for dd in q: disks[dd] = dd
        for i in range(4):
            nq = q[:i] + (swap(q, i),) + q[i+1:]
            k = tuple(sorted(nq))
            if k in seen: continue
            seen.add(k); nxt.append(nq)
    frontier = nxt
print("     ", len(quads), "Descartes configurations,", len(disks),
      "distinct disks, every one exact and null")
print()

def gsqrt(u2, r):
    """Exact sqrt in Z[i].  u=m+ni, u^2=(m^2-n^2)+2mn i, |u|^2=m^2+n^2=r.
    So m^2=(r+p)/2 and n^2=(r-p)/2 -- O(1), no search."""
    p, q = u2
    if (r + p) % 2 or (r - p) % 2: return None
    m2, n2 = (r + p)//2, (r - p)//2
    if m2 < 0 or n2 < 0: return None
    m, n = isqrt(m2), isqrt(n2)
    if m*m != m2 or n*n != n2: return None
    if q < 0: n = -n
    return (m, n) if gmul((m,n),(m,n)) == u2 else None

print("CLAIM 2  |spin(A,B)|^2 == curvature(A) + curvature(B)")
pairs, bad, sumsq = set(), 0, 0
for q in quads:
    for i in range(4):
        for j in range(i+1, 4):
            A, B = q[i], q[j]
            if (A, B) in pairs: continue
            pairs.add((A, B))
            u2 = gsub(gsmul(A[0], B[1]), gsmul(B[0], A[1]))
            n = gnorm(u2); r = isqrt(n)
            if r*r != n or r != A[0] + B[0]:
                bad += 1; continue
            u = gsqrt(u2, r)
            if u is None or u[0]**2 + u[1]**2 != A[0] + B[0]: bad += 1
            else: sumsq += 1
print("   tangent pairs tested :", len(pairs))
print("   failures             :", bad)
print("   pairs where b_A+b_B came out as m^2+n^2, m,n integers:", sumsq)
assert bad == 0
print("   => in an integral packing the sum of two touching curvatures is")
print("      ALWAYS a sum of two squares.  The spinor is the witness.")
print()
print("   sample:")
for (A, B) in sorted(pairs, key=lambda p: (p[0][0]+p[1][0], p[0][0]))[:9]:
    u2 = gsub(gsmul(A[0], B[1]), gsmul(B[0], A[1]))
    r = isqrt(gnorm(u2)); u = gsqrt(u2, r)
    print("      %5d + %5d = %6d   spinor %4d%+di   ->  %d^2 + %d^2"
          % (A[0], B[0], A[0]+B[0], u[0], u[1], u[0], u[1]))
