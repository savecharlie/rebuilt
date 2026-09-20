"""PAMGuard HyperbolicLocaliser.calcErrors, as written vs as intended.  Iris, fire 282.
As written: TOADInformation.clone() is super.clone() (shallow), so the jitter loop
writes into the ORIGINAL delay array, cumulatively; processTOADs3D is then called on
the original.  Net effect: each bootstrap sample sits on a random walk."""
import numpy as np, gs_check as g
rng = np.random.default_rng(282)
N, sig = 100, 1e-4           # bootStrapN default; a 0.1 ms configured TOAD error
ratios = []
for s in g.W[:40]:
    t0 = np.linalg.norm(g.G-s,axis=1)/g.C
    D0 = t0[None,:]-t0[:,None]                 # delays[i][j] = t_j - t_i
    def solve(D):                               # feed delays straight in (i<j used)
        return g.pamguard_gs(g.G, D[0])[0] if False else _solveD(D)
    def _solveD(D):
        # reconstruct arrival times relative to phone 0 from row 0, matching i<j use closely enough:
        # instead, solve with every pair using D directly
        n=len(g.G); cen=g.G.mean(0); L=[];R=[]
        for i in range(n):
            for j in range(i+1,n):
                d=-D[i][j]*g.C; row=np.zeros(3+n-1); row[:3]=g.G[i]-g.G[j]; row[3+i]=d
                rhs=0.5*(d*d+np.sum((g.G[i]-cen)**2)-np.sum((g.G[j]-cen)**2)); w=np.linalg.norm(g.G[i]-g.G[j])
                L.append(row*w); R.append(rhs*w)
        return cen+np.linalg.lstsq(np.array(L),np.array(R),rcond=None)[0][:3]
    # intended: independent jitter each sample
    P_int = np.array([_solveD(D0 + rng.normal(0,sig,D0.shape)) for _ in range(N)])
    # as written: shared array, cumulative
    D = D0.copy(); P_wr=[]
    for _ in range(N):
        D += rng.normal(0,sig,D.shape); P_wr.append(_solveD(D))
    P_wr=np.array(P_wr)
    ratios.append(np.linalg.norm(P_wr.std(0))/np.linalg.norm(P_int.std(0)))
    drift=np.linalg.norm(P_wr[-1]-s)
r=np.array(ratios)
print(f"reported-error inflation (as written / intended), 40 whales: median {np.median(r):.1f}x, range {r.min():.1f}-{r.max():.1f}x")
print(f"delay array after one call: drifted by ~sqrt(100)*sigma = {10*sig*1e3:.1f} ms per entry (was exact)")
