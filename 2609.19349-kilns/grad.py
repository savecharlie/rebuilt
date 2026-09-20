"""Does the gradiometer alone explain the hybrid data decaying faster than the
synthetic Bz?  Same sheet kiln as check.py, direct sum; vertical gradient by
central difference (dz=5 mm) and a 0.65 m two-sensor difference."""
import numpy as np
from check import bz_down, disk_cells, M, T
g = np.linspace(-4, 4, 81); X, Y = np.meshgrid(g, g, indexing="ij")
cx, cy = disk_cells(0.05); src = [(a, b, 0.8) for a, b in zip(cx, cy)]
ms = M * 0.05**2 * T
B = lambda h: bz_down((X, Y), h, src, ms)
HS = [0.2, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0]
bz, dg, gd = [], [], []
for h in HS:
    bz.append(B(h).max())
    dg.append(((B(h - .005) - B(h + .005)) / .01).max())
    gd.append((B(h) - B(h + .65)).max())
hyb1 = [100, 57.7, 27.2, 14.7, 8.8, 4.1, 2.4, 1.4]; hyb2 = [100, 59.7, 31.7, 19.3, 13.8, 7.4, 4.3, 2.7]
print(f"{'h':>4} {'Bz%':>6} {'dBz/dz%':>8} {'0.65m-diff%':>11} | {'hyb P1':>6} {'hyb P2':>6}")
for i, h in enumerate(HS):
    print(f"{h:4.1f} {100*bz[i]/bz[0]:6.1f} {100*dg[i]/dg[0]:8.1f} {100*gd[i]/gd[0]:11.1f} | {hyb1[i]:6.1f} {hyb2[i]:6.1f}")
