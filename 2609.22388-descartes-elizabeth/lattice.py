"""Which Gaussian integers occur as tangency spinors of the (-1,2,2,3) packing?
And is there a constraint on which curvature-sums can occur at all?"""
from math import isqrt, gcd

def gmul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def gadd(a,b): return (a[0]+b[0], a[1]+b[1])
def gsub(a,b): return (a[0]-b[0], a[1]-b[1])
def gsmul(k,a): return (k*a[0], k*a[1])
def gnorm(a): return a[0]*a[0]+a[1]*a[1]

def swap(q,i):
    sb = sum(d[0] for d in q) - q[i][0]
    sw = (0,0)
    for j,d in enumerate(q):
        if j != i: sw = gadd(sw, d[1])
    return (2*sb - q[i][0], gsub(gsmul(2,sw), q[i][1]))

def gsqrt(u2, r):
    p,q = u2
    if (r+p)%2 or (r-p)%2: return None
    m2,n2 = (r+p)//2,(r-p)//2
    if m2<0 or n2<0: return None
    m,n = isqrt(m2),isqrt(n2)
    if m*m!=m2 or n*n!=n2: return None
    if q<0: n=-n
    return (m,n) if gmul((m,n),(m,n))==u2 else None

seed = ((-1,(0,0)), (2,(-1,0)), (2,(1,0)), (3,(0,2)))
seen = {tuple(sorted(seed))}; quads, frontier = [], [seed]
for _ in range(7):
    nxt = []
    for q in frontier:
        quads.append(q)
        for i in range(4):
            nq = q[:i] + (swap(q,i),) + q[i+1:]
            k = tuple(sorted(nq))
            if k in seen: continue
            seen.add(k); nxt.append(nq)
    frontier = nxt

spins, pairs, disks = set(), set(), set()
for q in quads:
    for d in q: disks.add(d)
    for i in range(4):
        for j in range(i+1,4):
            A,B = q[i],q[j]
            if (A,B) in pairs: continue
            pairs.add((A,B))
            u2 = gsub(gsmul(A[0],B[1]), gsmul(B[0],A[1]))
            r = isqrt(gnorm(u2)); u = gsqrt(u2,r)
            assert u and u[0]**2+u[1]**2 == A[0]+B[0], (A,B,u)
            m,n = abs(u[0]),abs(u[1])
            spins.add((max(m,n),min(m,n)))
print("configurations:", len(quads), " disks:", len(disks), " tangent pairs:", len(pairs))
print("distinct spinors up to units and conjugation:", len(spins))

nc = sorted(s for s in spins if gcd(s[0],s[1]) != 1)
print("spinors with gcd(m,n) != 1:", len(nc), nc[:12])

def factors(n):
    n=abs(n); out=[]; d=2
    while d*d<=n:
        while n%d==0: out.append(d); n//=d
        d+=1
    if n>1: out.append(n)
    return out

sums = sorted({A[0]+B[0] for (A,B) in pairs})
bad = {s:[p for p in factors(s) if p%4==3] for s in sums}
bad = {k:v for k,v in bad.items() if v}
print("distinct curvature-sums:", len(sums), "range", sums[0], "..", sums[-1])
print("with a prime factor = 3 (mod 4):", len(bad), list(bad.items())[:5])
print()
print("smallest 26 curvature-sums:", sums[:26])
print("missing small integers (not a curvature-sum):",
      [n for n in range(1, 60) if n not in set(sums)])
