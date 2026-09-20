"""Rebuild Hegyi et al. 2609.19349 Table 1 (synthetic kiln) and check it against
direct dipole summation with no FFT.  Iris the Maker, fire 285.

Their recipe (Sec 3.2.1-3.2.4):
  grid 101x101, dx=0.1 m, 10 m square; cells with centres inside r<=1.0 m get
  m_sub = M dx^2 T = 1.5 * 0.01 * 1.0 = 0.015 A m^2, aligned I=70, D=3,
  ALL at one depth d_sub = 0.8 m.  Bz computed on z=0, then +1 nT/m trend in x,
  heterogeneity (0.5 nT, smoothed 10 px), noise, then FFT upward continuation
  exp(-|k| h).  Peak = global max.  Retention relative to h=0.2.
"""
import numpy as np
from scipy.ndimage import gaussian_filter

MU = 1e-7
M, T, R = 1.5, 1.0, 1.0
I, Dec = np.radians(70), np.radians(3)
# unit vector of magnetisation, geophysics frame x=east? use x=north,y=east,z=down
u = np.array([np.cos(I) * np.cos(Dec), np.cos(I) * np.sin(Dec), np.sin(I)])


def bz_down(obs_xy, h, src_xyz, msub):
    """vertical (down-positive) field in nT at height h above ground, from point
    dipoles at src_xyz (z = depth, down positive)."""
    X, Y = obs_xy
    out = np.zeros_like(X)
    for sx, sy, sz in src_xyz:
        rx, ry, rz = X - sx, Y - sy, (-h) - sz      # z down: obs at -h
        r2 = rx * rx + ry * ry + rz * rz
        r = np.sqrt(r2)
        mdotr = msub * (u[0] * rx + u[1] * ry + u[2] * rz)
        out += MU * (3 * mdotr * rz / r ** 5 - msub * u[2] / r ** 3)
    return out * 1e9


def disk_cells(dx=0.1):
    g = np.arange(-5, 5 + 1e-9, dx)
    cx, cy = np.meshgrid(g, g, indexing="ij")
    mask = cx ** 2 + cy ** 2 <= R ** 2
    return cx[mask], cy[mask]


def their_pipeline(extras=True, seed=0):
    g = np.linspace(-5, 5, 101)
    X, Y = np.meshgrid(g, g, indexing="ij")
    cx, cy = disk_cells()
    src = [(a, b, 0.8) for a, b in zip(cx, cy)]
    B0 = bz_down((X, Y), 0.0, src, M * 0.01 * T)
    if extras:
        rng = np.random.default_rng(seed)
        B0 = B0 + 1.0 * X                                  # trend along x
        B0 = B0 + gaussian_filter(rng.normal(0, 0.5, X.shape), 10)
        B0 = B0 + gaussian_filter(rng.normal(0, 0.00894, X.shape), 5)
    k = 2 * np.pi * np.fft.fftfreq(101, d=0.1)
    KX, KY = np.meshgrid(k, k, indexing="ij")
    K = np.hypot(KX, KY)
    F = np.fft.fft2(B0)
    return {h: np.real(np.fft.ifft2(F * np.exp(-K * h))) for h in HS}, B0


def truth(volume=False, n=6):
    """direct summation on a 60 m square, no periodicity. volume=True spreads the
    dipoles through the 1 m thickness (top 0.3 m) instead of one sheet at 0.8 m."""
    g = np.linspace(-6, 6, 121)          # peak is near centre; 12 m is plenty
    X, Y = np.meshgrid(g, g, indexing="ij")
    cx, cy = disk_cells(0.05)
    if volume:
        zs = 0.3 + (np.arange(n) + 0.5) * T / n
        src = [(a, b, z) for z in zs for a, b in zip(cx, cy)]
        ms = M * 0.05 ** 2 * T / n
    else:
        src = [(a, b, 0.8) for a, b in zip(cx, cy)]
        ms = M * 0.05 ** 2 * T
    return {h: bz_down((X, Y), h, src, ms) for h in HS}


HS = [0.2, 0.3, 0.4, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 10.0, 15.0]
PAPER = [322.7, 277.5, 239.0, 206.4, 104.5, 58.4, 35.8, 23.7, 16.8, 12.7, 10.1, 8.4, 7.3, 6.0, 4.8, 4.8]

if __name__ == "__main__":
    fft_x, B0 = their_pipeline(True)
    fft_c, B0c = their_pipeline(False)
    tru = truth(False)
    vol = truth(True)
    print(f"window mean of clean Bz(z=0) on their 10 m grid: {B0c.mean():.2f} nT")
    print(f"window mean with trend+noise:                    {B0.mean():.2f} nT")
    print(f"{'h':>5} {'paper':>7} {'fft+ex':>7} {'fftcln':>7} {'direct':>7} {'volume':>7} | "
          f"{'ret_p':>6} {'ret_d':>6} {'ret_v':>6}")
    rows = []
    for h, p in zip(HS, PAPER):
        a, b, c, d = fft_x[h].max(), fft_c[h].max(), tru[h].max(), vol[h].max()
        rows.append((h, p, a, b, c, d))
    p0, c0, d0 = rows[0][1], rows[0][4], rows[0][5]
    for h, p, a, b, c, d in rows:
        print(f"{h:5.1f} {p:7.1f} {a:7.1f} {b:7.1f} {c:7.2f} {d:7.2f} | "
              f"{100*p/p0:6.1f} {100*c/c0:6.2f} {100*d/d0:6.2f}")
    np.save("rows.npy", np.array(rows))
