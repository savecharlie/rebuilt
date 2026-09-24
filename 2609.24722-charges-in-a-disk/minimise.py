"""Road 3: my own minimum-energy configurations, so the fit is not borrowed.

Nothing here is trusted until it reproduces energies somebody else published:
  N = 2, 3, 4   closed forms
  N = 60, 61, 92, 99   Lavrov & Nikonov 2026 (arXiv:2609.20777) global-minimum
                       candidates, obtained by a different method (quenched MD)

Parameterisation follows the paper: Nb charges pinned to the rim with free
angle, the rest interior with r = sigma sin^2(t), so the hard wall is exact
and no penalty term is needed.
"""
import numpy as np
from scipy.optimize import minimize

RNG = np.random.default_rng


def nb_formula(N):
    return 2.84328 * N ** (2 / 3) - 0.530196 * N ** (1 / 3) - 2.32866


def cart_energy_grad(xy):
    """E = sum_{i<j} 1/r_ij and dE/d(x,y), via BLAS.

    r2 = |xi|^2 + |xj|^2 - 2 xi.xj is one dgemm, and then both the row sums
    and the force accumulation are matrix products too:
        g_i = (W x)_i - (sum_j W_ij) x_i,   W_ij = r_ij^{-3}
    The Gram form loses digits when two points are close, so `slow_energy`
    below is kept as the control and they are compared in validate().
    """
    n = xy.shape[0]
    sq = np.einsum("ij,ij->i", xy, xy)
    r2 = sq[:, None] + sq[None, :] - 2.0 * (xy @ xy.T)
    np.fill_diagonal(r2, np.inf)
    np.maximum(r2, 1e-300, out=r2)
    inv = r2 ** -0.5
    E = 0.5 * inv.sum()
    W = inv ** 3
    s = W.sum(axis=1)
    g = W @ xy - s[:, None] * xy
    return E, g


def slow_energy(xy, block=512):
    """direct pairwise energy, no cancellation — the control for the above"""
    n = xy.shape[0]
    E = 0.0
    for a in range(0, n, block):
        b = min(a + block, n)
        d = xy[a:b, None, :] - xy[None, :, :]
        r2 = np.einsum("ijk,ijk->ij", d, d)
        r2[np.arange(b - a), np.arange(a, b)] = np.inf
        E += (r2 ** -0.5).sum()
    return 0.5 * E


def pack(t, ui, ub):
    return np.concatenate([t, ui, ub])


def unpack(v, ni):
    return v[:ni], v[ni:2 * ni], v[2 * ni:]


def positions(t, ui, ub, sigma):
    r = sigma * np.sin(t) ** 2
    xi = np.stack([r * np.cos(ui), r * np.sin(ui)], axis=1)
    xb = np.stack([np.cos(ub), np.sin(ub)], axis=1)
    return np.vstack([xi, xb]), r


def objective(v, ni, sigma):
    t, ui, ub = unpack(v, ni)
    xy, r = positions(t, ui, ub, sigma)
    E, g = cart_energy_grad(xy)
    gi, gb = g[:ni], g[ni:]
    dt = (gi[:, 0] * np.cos(ui) + gi[:, 1] * np.sin(ui)) * sigma * np.sin(2 * t)
    dui = -gi[:, 0] * r * np.sin(ui) + gi[:, 1] * r * np.cos(ui)
    dub = -gb[:, 0] * np.sin(ub) + gb[:, 1] * np.cos(ub)
    return E, np.concatenate([dt, dui, dub])


def sample_equilibrium(n, rng):
    """r with density proportional to the conducting-disk measure: r = sin(t),
    t uniform-by-weight sin t  ->  inverse cdf r = sqrt(1-(1-U)^2)? no:
    P(R<r) = 1 - sqrt(1-r^2), so r = sqrt(1-(1-U)^2)."""
    u = rng.random(n)
    return np.sqrt(1.0 - (1.0 - u) ** 2)


def relax(N, Nb, seed=0, ramp=(0.0, 0.5, 0.8, 1.0), maxiter=4000):
    rng = RNG(seed)
    ni = N - Nb
    sigma0 = 1.0 - 1.0 / (2.0 * np.sqrt(N))
    r0 = np.clip(sample_equilibrium(ni, rng) * sigma0, 1e-4, sigma0 - 1e-6)
    t = np.arcsin(np.sqrt(np.clip(r0 / sigma0, 0, 1)))
    ui = rng.random(ni) * 2 * np.pi
    ub = np.sort(rng.random(Nb)) * 2 * np.pi if Nb else np.zeros(0)
    if Nb:
        ub = np.linspace(0, 2 * np.pi, Nb, endpoint=False) + rng.random() * 0.1
    v = pack(t, ui, ub)
    for frac in ramp:
        sigma = sigma0 + frac * (1.0 - sigma0)
        res = minimize(objective, v, args=(ni, sigma), jac=True,
                       method="L-BFGS-B",
                       options=dict(maxiter=maxiter, maxfun=maxiter * 2,
                                    ftol=1e-16, gtol=1e-12))
        v = res.x
    t, ui, ub = unpack(v, ni)
    xy, _ = positions(t, ui, ub, 1.0)
    E, _ = cart_energy_grad(xy)
    return E, xy


def best_of(N, seeds=8, nb_scan=(0,), maxiter=4000):
    best, bxy, bnb = np.inf, None, None
    nb0 = int(round(nb_formula(N)))
    for dnb in nb_scan:
        Nb = max(0, min(N, nb0 + dnb))
        for s in range(seeds):
            E, xy = relax(N, Nb, seed=1000 * (dnb + 5) + s, maxiter=maxiter)
            if E < best:
                best, bxy, bnb = E, xy, Nb
    return best, bxy, bnb


PUBLISHED = {60: 2159.3584240930, 61: 2237.19264190,
             92: 5358.35353314, 99: 6254.83029083}
CLOSED = {2: 0.5, 3: np.sqrt(3.0), 4: 2 * np.sqrt(2.0) + 1.0}


def validate():
    print("=" * 74)
    print("CONTROL — fast BLAS kernel vs direct pairwise sum, and the gradient")
    print("=" * 74)
    rng = RNG(7)
    for n in (200, 1000, 2500):
        xy = rng.random((n, 2)) * 1.6 - 0.8
        E, g = cart_energy_grad(xy)
        Es = slow_energy(xy)
        # finite difference on one coordinate
        h, k = 1e-6, 5
        xp = xy.copy(); xp[k, 0] += h
        xm = xy.copy(); xm[k, 0] -= h
        fd = (slow_energy(xp) - slow_energy(xm)) / (2 * h)
        print(f"  n={n:5d}  E rel diff {abs(E-Es)/Es:.2e}   "
              f"grad[{k},0] {g[k,0]:+.8f} vs fd {fd:+.8f}  "
              f"rel {abs(g[k,0]-fd)/abs(fd):.2e}")
    print()
    print("=" * 74)
    print("CONTROL — closed forms")
    print("=" * 74)
    for N, want in CLOSED.items():
        E, _, nb = best_of(N, seeds=4, nb_scan=(0,), maxiter=2000)
        Eb, _, _ = best_of(N, seeds=4, nb_scan=(N - int(round(nb_formula(N))),))
        E = min(E, Eb)
        print(f"  N={N:3d}  mine {E:.10f}   exact {want:.10f}   "
              f"rel {abs(E-want)/want:.2e}")
    print()
    print("=" * 74)
    print("CONTROL — published global-minimum candidates (arXiv:2609.20777)")
    print("=" * 74)
    for N, want in PUBLISHED.items():
        E, _, nb = best_of(N, seeds=14, nb_scan=(-1, 0, 1))
        print(f"  N={N:3d}  mine {E:.8f}   published {want:.8f}   "
              f"rel {(E-want)/want:+.2e}   Nb={nb}")


if __name__ == "__main__":
    validate()
