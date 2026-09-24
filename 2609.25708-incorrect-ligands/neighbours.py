"""Second exact pass: the neighbour-resolved ensemble.

Their Table I gives r_off^(a) = r_off^(0) exp(-d_eps * a), and the binding-state
coupling shifts a receptor's OWN binding energy by -4 J_b (k-2), where k is its
number of BOUND nearest neighbours (a pinned boundary neighbour counts as free,
since it is the absence of a receptor).  So

    K_d(k, a) = r_off/r_on = exp( eps_b + d_eps_b - 4 J_b (k-2) - d_eps * a )

and at k = 2 this returns exactly exp(-1) = 0.3679, which is Table I's quoted
r_off^(0) = 0.37 /t0.  That is an internal check on the whole reading.

A bound receptor reaches full activation only if it survives TWO enzyme steps at
rate r_a against unbinding -- their kinetic proofreading:

    P_act(k) = [ r_a / (r_a + K_d(k,0)) ] * [ r_a / (r_a + K_d(k,1/2)) ]

So the expected signal of the seed-free detector is  <sum_i bound_i P_act(k_i)>,
which is an EQUILIBRIUM average I can take exactly.

The reweighting factor depends only on (B, P, n, s_c), so conditional sums of the
neighbour-count populations N_k inside each cell are fixed -- store them once.
"""
import numpy as np, os, math
from ising import geometry

def build(L=5, chunk_bits=20, cache="nbr_L%d.npz"):
    path = cache % L
    if os.path.exists(path):
        z = np.load(path); return z["H"], z["S"], (int(z["Bmin"]), int(z["Pmin"]))
    bonds, miss = geometry(L); n_sites = L*L
    nbr = [[] for _ in range(n_sites)]
    for i, j in bonds: nbr[i].append(j); nbr[j].append(i)
    Bmin, Pmin = -len(bonds), -int(miss.sum())
    nB, nP, nN = 2*len(bonds)+1, 2*int(miss.sum())+1, n_sites+1
    centre = (L//2)*L + (L//2)
    H = np.zeros(nB*nP*nN*2, dtype=np.int64)
    S = np.zeros((5, nB*nP*nN*2), dtype=np.int64)     # S[k] = sum of N_k over cell
    total = 1 << n_sites; step = 1 << min(chunk_bits, n_sites)
    for start in range(0, total, step):
        idx = np.arange(start, min(start+step, total), dtype=np.int64)
        up = np.empty((n_sites, idx.size), dtype=np.int8)
        for q in range(n_sites): up[q] = ((idx >> q) & 1).astype(np.int8)
        s = up.astype(np.int32)*2 - 1
        B = np.zeros(idx.size, np.int32)
        for i, j in bonds: B += s[i]*s[j]
        P = np.zeros(idx.size, np.int32)
        for i in range(n_sites):
            if miss[i]: P += int(miss[i])*s[i]
        nUp = up.sum(axis=0).astype(np.int32)
        key = (((B-Bmin).astype(np.int64)*nP + (P-Pmin))*nN + nUp)*2 + up[centre]
        H += np.bincount(key, minlength=H.size)
        kcount = np.zeros((n_sites, idx.size), np.int8)
        for i in range(n_sites):
            for j in nbr[i]: kcount[i] += up[j]
        for k in range(5):
            tally = np.zeros(idx.size, np.int32)
            for i in range(n_sites):
                tally += (up[i] == 1) & (kcount[i] == k)
            S[k] += np.bincount(key, weights=tally, minlength=H.size).astype(np.int64)
    sh = (nB, nP, nN, 2)
    H = H.reshape(sh); S = S.reshape((5,)+sh)
    assert H.sum() == total
    np.savez_compressed(path, H=H, S=S, Bmin=Bmin, Pmin=Pmin)
    return H, S, (Bmin, Pmin)


class Detector:
    EPS_B, D_EPS_B, D_EPS = -3.5, 2.5, 4.0
    R_A = 0.042
    def __init__(self, L=5):
        self.L = L; self.H, self.S, (self.Bmin, self.Pmin) = build(L)
        nB, nP, nN, _ = self.H.shape
        self.Bv = np.arange(self.Bmin, self.Bmin+nB)[:, None, None, None]
        self.Pv = np.arange(self.Pmin, self.Pmin+nP)[None, :, None, None]
        self.nv = np.arange(nN)[None, None, :, None]
        self.scv = np.array([-1, 1])[None, None, None, :]
        self.M = 2*self.nv - L*L
        self.nz = self.H > 0
    def _w(self, J, h, h_seed=0.0):
        E = -J*self.Bv + J*self.Pv - h*self.M - h_seed*self.scv
        lw = np.where(self.nz, np.log(np.where(self.nz, self.H, 1)) - E, -np.inf)
        return np.exp(lw - lw.max())
    def Nk(self, J, h, h_seed=0.0):
        w = self._w(J, h, h_seed); Z = w.sum()
        return np.array([(w*self.S[k]/np.where(self.nz, self.H, 1)).sum()/Z for k in range(5)])
    def mean_n(self, J, h, h_seed=0.0):
        w = self._w(J, h, h_seed); return float((w*self.nv).sum()/w.sum())
    def Kd(self, J, k, a):
        return math.exp(self.EPS_B + self.D_EPS_B - 4*J*(k-2) - self.D_EPS*a)
    def P_act(self, J, k):
        r = self.R_A
        return (r/(r+self.Kd(J,k,0.0))) * (r/(r+self.Kd(J,k,0.5)))
    def signal(self, J, h, h_seed=0.0):
        """<number of fully activated receptors> in the quasi-static neighbourhood
        approximation.  SEAM: treats each bound receptor's k as frozen over the
        two enzyme steps."""
        Nk = self.Nk(J, h, h_seed)
        return float(sum(Nk[k]*self.P_act(J, k) for k in range(5))), Nk
