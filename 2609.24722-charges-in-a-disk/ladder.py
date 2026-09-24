"""Generate my own E(N) ladder, seeded by a warped triangular lattice.

The seed matters more than the optimiser here: a triangular lattice pushed
through the map that carries the uniform measure to the equilibrium measure,

        r = sqrt(1 - (1 - r_ref^2)^2),

starts within a fraction of a per cent of the minimum, so L-BFGS only has to
polish.  Random starts need thousands of iterations to do the same job.
"""
import json, sys, time
import numpy as np
from scipy.optimize import minimize
from minimise import (cart_energy_grad, slow_energy, nb_formula, objective,
                      pack, unpack, positions, RNG)


def warped_lattice(ni, rng, sigma):
    """ni points, triangular lattice warped to the equilibrium density"""
    a = np.sqrt(2 * np.pi / (np.sqrt(3.0) * max(ni, 1)))
    k = int(np.ceil(1.6 / a)) + 2
    i, j = np.meshgrid(np.arange(-k, k + 1), np.arange(-k, k + 1), indexing="ij")
    x = a * (i + 0.5 * j).ravel()
    y = a * (np.sqrt(3.0) / 2.0 * j).ravel()
    th = rng.random() * np.pi / 3
    xr = x * np.cos(th) - y * np.sin(th) + (rng.random() - .5) * a
    yr = x * np.sin(th) + y * np.cos(th) + (rng.random() - .5) * a
    r = np.hypot(xr, yr)
    keep = np.argsort(r)[:ni]
    xr, yr, r = xr[keep], yr[keep], r[keep]
    r = np.clip(r / max(r.max(), 1e-12), 0, 1)
    rw = np.sqrt(np.clip(1.0 - (1.0 - r ** 2) ** 2, 0, 1)) * sigma
    ang = np.arctan2(yr, xr)
    return np.clip(rw, 1e-6, sigma * (1 - 1e-9)), ang


def run(N, Nb, seed, maxiter, ramp=(0.0, 1.0)):
    rng = RNG(seed)
    ni = N - Nb
    sigma0 = 1.0 - 1.0 / (2.0 * np.sqrt(N))
    r0, ui = warped_lattice(ni, rng, sigma0)
    t = np.arcsin(np.sqrt(np.clip(r0 / sigma0, 0, 1)))
    ub = np.linspace(0, 2 * np.pi, Nb, endpoint=False) + rng.random() * 0.3
    v = pack(t, ui, ub)
    for frac in ramp:
        sigma = sigma0 + frac * (1.0 - sigma0)
        res = minimize(objective, v, args=(ni, sigma), jac=True,
                       method="L-BFGS-B",
                       options=dict(maxiter=maxiter, maxfun=2 * maxiter,
                                    ftol=1e-17, gtol=1e-13))
        v = res.x
    t, ui, ub = unpack(v, ni)
    xy, _ = positions(t, ui, ub, 1.0)
    E, _ = cart_energy_grad(xy)
    return E, xy


def best(N, seeds, maxiter, dnb=(-1, 0, 1)):
    nb0 = int(round(nb_formula(N)))
    bE, bxy, bNb = np.inf, None, None
    for d in dnb:
        Nb = max(3, min(N - 1, nb0 + d))
        for s in range(seeds):
            E, xy = run(N, Nb, 97 * s + 13 * (d + 3), maxiter)
            if E < bE:
                bE, bxy, bNb = E, xy, Nb
    return bE, bxy, bNb


def main(path="ladder.json"):
    # (N, seeds, maxiter).  The warped-lattice seed starts close, so a
    # couple of restarts is enough; the Nb scan matters more than the seed.
    plan = [(100, 3, 2500), (150, 3, 2500), (220, 2, 2000), (320, 2, 2000),
            (460, 2, 1600), (660, 2, 1400), (950, 1, 1200), (1350, 1, 1000),
            (1900, 1, 800)]
    out = []
    for N, seeds, mx in plan:
        t0 = time.time()
        E, xy, Nb = best(N, seeds, mx)
        Es = slow_energy(xy)                      # control: no Gram cancellation
        r = np.hypot(xy[:, 0], xy[:, 1])
        rec = dict(N=N, Nb=Nb, E=E, E_direct=Es,
                   rel_kernel=abs(E - Es) / Es,
                   rmax=float(r.max()), secs=round(time.time() - t0, 1))
        out.append(rec)
        print(f"N={N:5d} Nb={Nb:4d}  E={E:.10f}  kernel rel "
              f"{rec['rel_kernel']:.1e}  rmax={rec['rmax']:.6f}  "
              f"{rec['secs']:6.1f}s", flush=True)
        json.dump(out, open(path, "w"), indent=1)
    print("done", flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:])
