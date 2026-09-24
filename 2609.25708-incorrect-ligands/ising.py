"""
Exact equilibrium statistical mechanics of the receptor cluster of
arXiv:2609.25708v1 (Chen, Leung & Chen), using the authors' OWN mapping.

Their eq. (1) with all activation levels a_i = 0, under s_i = 2 b_i^2 - 1:

    H = -J sum_<ij> s_i s_j - sum_i h_i s_i ,   h = (ln c - eps_b - d_eps_b)/2

eps_b = -3.5, d_eps_b = 2.5  =>  h = (ln c + 1)/2,  h = 0 at c* = e^-1.

End Matter: a boundary receptor's absent outside neighbour is held at b = 0, i.e.
s = -1.  A pinned neighbour at s_p adds -J s_i s_p to H, i.e. a LOCAL FIELD of
+J*s_p = -J.  So site i carries an extra field -J * m_i, m_i = missing neighbours.

KEY TRICK.  With h_i = h + (-J m_i) + h_seed*delta_{i,centre},

    E = -J*B  +  J*P  -  h*M  -  h_seed*s_c

where B = sum_<ij> s_i s_j, P = sum_i m_i s_i, M = sum_i s_i, s_c = centre spin.
Four integers.  So ONE pass over all 2^25 states, histogrammed by (B, P, n, s_c),
answers every (J, h, h_seed) afterwards by reweighting.  Nothing is re-simulated.

Iris (Opus 5), fire 298.
"""
import numpy as np, itertools, os, json

L = 5
N = L * L
CENTRE = (L // 2) * L + (L // 2)          # row-major index of the middle receptor


def geometry(L):
    """bonds, missing-neighbour count m_i, for an LxL free cluster."""
    bonds = []
    for r in range(L):
        for cc in range(L):
            i = r * L + cc
            if cc + 1 < L: bonds.append((i, i + 1))
            if r + 1 < L:  bonds.append((i, i + L))
    m = np.zeros(L * L, dtype=np.int64)
    for r in range(L):
        for cc in range(L):
            i = r * L + cc
            m[i] = (cc == 0) + (cc == L - 1) + (r == 0) + (r == L - 1)
    return bonds, m


BONDS, MISS = geometry(L)
assert len(BONDS) == 40 and MISS.sum() == 20


# ---------------------------------------------------------------- enumeration
def enumerate_histogram(L=5, chunk_bits=20, cache="hist_L%d.npz"):
    """Exact histogram over all 2^(L*L) states, keyed by (B, P, n, s_c)."""
    path = cache % L
    if os.path.exists(path):
        z = np.load(path)
        return z["H"], tuple(z["axes"])
    bonds, miss = geometry(L)
    n_sites = L * L
    Bmin, Bmax = -len(bonds), len(bonds)
    Pmin, Pmax = -int(miss.sum()), int(miss.sum())
    nB, nP, nN = Bmax - Bmin + 1, Pmax - Pmin + 1, n_sites + 1
    centre = (L // 2) * L + (L // 2)
    H = np.zeros(nB * nP * nN * 2, dtype=np.int64)

    total = 1 << n_sites
    step = 1 << min(chunk_bits, n_sites)
    for start in range(0, total, step):
        idx = np.arange(start, min(start + step, total), dtype=np.int64)
        s = np.empty((n_sites, idx.size), dtype=np.int8)
        for k in range(n_sites):
            s[k] = (((idx >> k) & 1) * 2 - 1).astype(np.int8)
        B = np.zeros(idx.size, dtype=np.int32)
        for (i, j) in bonds:
            B += s[i].astype(np.int32) * s[j]
        P = np.zeros(idx.size, dtype=np.int32)
        for i in range(n_sites):
            if miss[i]:
                P += int(miss[i]) * s[i].astype(np.int32)
        nUp = np.zeros(idx.size, dtype=np.int32)
        for k in range(n_sites):
            nUp += (s[k] > 0)
        sc = (s[centre] > 0).astype(np.int64)
        key = (((B - Bmin).astype(np.int64) * nP + (P - Pmin)) * nN + nUp) * 2 + sc
        H += np.bincount(key, minlength=H.size)
    H = H.reshape(nB, nP, nN, 2)
    assert H.sum() == total, (H.sum(), total)
    np.savez_compressed(path, H=H, axes=np.array([Bmin, Pmin, 0]))
    return H, (Bmin, Pmin, 0)


class Cluster:
    """Reweight the exact histogram to any (J, h, h_seed)."""

    def __init__(self, L=5):
        self.L = L
        self.H, (self.Bmin, self.Pmin, _) = enumerate_histogram(L)
        nB, nP, nN, _ = self.H.shape
        self.Bv = np.arange(self.Bmin, self.Bmin + nB)[:, None, None, None]
        self.Pv = np.arange(self.Pmin, self.Pmin + nP)[None, :, None, None]
        self.nv = np.arange(nN)[None, None, :, None]
        self.scv = np.array([-1, 1])[None, None, None, :]
        self.M = 2 * self.nv - L * L
        self.nz = self.H > 0

    def logw(self, J, h, h_seed=0.0):
        """log Boltzmann weight of each (B,P,n,s_c) cell, +log multiplicity."""
        E = -J * self.Bv + J * self.Pv - h * self.M - h_seed * self.scv
        lw = np.where(self.nz, np.log(np.where(self.nz, self.H, 1)) - E, -np.inf)
        return lw

    def _agg(self, J, h, h_seed=0.0):
        lw = self.logw(J, h, h_seed)
        mx = lw.max()
        w = np.exp(lw - mx)
        return w, mx

    def logZ(self, J, h, h_seed=0.0):
        w, mx = self._agg(J, h, h_seed)
        return np.log(w.sum()) + mx

    def mean_bound(self, J, h, h_seed=0.0):
        """<n_bound> / N."""
        w, _ = self._agg(J, h, h_seed)
        return float((w * self.nv).sum() / w.sum()) / (self.L * self.L)

    def Fn(self, J, h, h_seed=0.0):
        """Projected free energy F(n) = -ln sum_{states with n bound} e^-E, shifted to min 0."""
        w, mx = self._agg(J, h, h_seed)
        Zn = w.sum(axis=(0, 1, 3))
        F = -(np.log(np.where(Zn > 0, Zn, np.nan)) + mx)
        return F - np.nanmin(F)

    def pin_centre_up(self):
        """Return a Cluster-like view with the centre spin forced to +1 (strong-binding seed)."""
        c = object.__new__(Cluster)
        c.__dict__.update(self.__dict__)
        c.H = self.H.copy()
        c.H[:, :, :, 0] = 0                 # drop every state with s_centre = -1
        c.nz = c.H > 0
        return c


# ------------------------------------------------- independent route: transfer matrix
def transfer_logZ(L, J, hsite):
    """Exact log Z by column-to-column transfer.  hsite: length L*L array of local fields.
    Independent of the enumeration code path (law 6)."""
    states = np.array(list(itertools.product([-1, 1], repeat=L)), dtype=np.int64)  # (2^L, L)
    intra = J * (states[:, :-1] * states[:, 1:]).sum(axis=1)                       # within-column bonds
    T = J * (states @ states.T)                                                    # between-column bonds
    v = None
    for col in range(L):
        hcol = hsite[col * L:(col + 1) * L] if False else hsite[np.arange(L) * L + col]
        w = intra + states @ hcol
        if v is None:
            v = w.copy()
        else:
            mx = v.max()
            v = np.log(np.exp(v - mx) @ np.exp(T)) + mx + w
    mx = v.max()
    return float(np.log(np.exp(v - mx).sum()) + mx)


def site_fields(L, J, h, h_seed=0.0):
    _, miss = geometry(L)
    hs = h - J * miss.astype(float)
    hs[(L // 2) * L + (L // 2)] += h_seed
    return hs


# ------------------------------------------------- Onsager (external known answer)
def onsager_logZ_per_site(K, n=4000):
    """-beta f for the square-lattice Ising model, exact (Onsager 1944)."""
    th = (np.arange(n) + 0.5) * np.pi / n
    c2, s2 = np.cosh(2 * K), np.sinh(2 * K)
    A = c2 ** 2 - s2 * (np.cos(th)[:, None] + np.cos(th)[None, :])
    return np.log(2.0) + 0.5 * np.log(A).mean()


def strip_logZ_per_site(W, Lc, K):
    """Transfer matrix for a W-wide strip, PERIODIC across the width, Lc columns,
    periodic along the length.  Used only to approach the Onsager limit."""
    states = np.array(list(itertools.product([-1, 1], repeat=W)), dtype=np.int64)
    intra = K * (states * np.roll(states, -1, axis=1)).sum(axis=1)
    # Tm is SYMMETRIC by construction; use the symmetric solver (exact, and far faster)
    Tm = np.exp(K * (states @ states.T) + 0.5 * (intra[:, None] + intra[None, :]))
    assert np.allclose(Tm, Tm.T, rtol=0, atol=0), "transfer matrix must be symmetric"
    lam = np.linalg.eigvalsh(Tm)
    return float(np.log(lam[-1]) / W), float(np.log(lam[-1] / lam[-2]))
