"""Tao, Esteves, Lee et al., arXiv 2609.15920 -- can the holographic heights be the
heights the cells are advected at?  Iris the Maker, fire 281.

The paper writes Vx = Vx_swim + Vx_flow with Vx_flow ~ -gamma*<z>, and says the Fig 3D
slopes 'are a measure of the average height ... in agreement with the results from
holography'.  Since a cell cannot swim upstream faster than it swims at all,
<Vx_swim> <= <Vswim>, which gives a hard ceiling on the advection height:

    <z_eff> <= (<Vswim> - <Vx>) / gamma

and it is sharpest at the highest shear.  Reached with EVERY cell pointed dead upstream;
the real value is lower (Fig 2E: 32% of WT move upstream at all).

Sources, all read off the paper by eye (Fig 3D at 400 dpi, Fig 2D, Fig 1C, Fig 2C):
heights  = Fig 3C means (holography, flat surface, window 0-5 um), SD and N_traj
slopes   = printed on Fig 3D
Vx(21.7) = Fig 2D normalized points x Vswim, Fig 3D points for the other species
Vswim    = Fig 1C (E. coli 19.8, P.a. 26.2, V.c. 64.2), Fig 2C read by eye (F ~21.2, A ~16.5)
S. enterica Vswim is in the SI, which is not on arXiv -> no bound for it.
"""
import math
G = 21.7
rows = [  # name, z_mean, z_sd, n_traj, printed_slope, Vx_at_G, Vswim
    ("S. enterica",    1.61, 1.22,  5, -0.62, -7,  None),
    ("E. coli PhiF",   1.54, 1.10, 33, None,  -0.44*21.2, 21.2),
    ("E. coli WT",     1.64, 0.77, 23, -0.57, -0.47*19.8, 19.8),
    ("E. coli PhiA",   2.34, 1.09, 25, None,  -1.22*16.5, 16.5),
    ("V. cholerae",    2.25, 1.49, 21, -2.42, -44, 64.2),
    ("P. aeruginosa",  2.74, 1.42,  8, -3.61, -70, 26.2),
]
print(f"{'strain':15s} {'<z> holo':>12s} {'|slope|':>8s} {'ceiling':>8s}  verdict")
for n, z, sd, N, s, vx, vs in rows:
    sem = sd / math.sqrt(N)
    ceil = (vs - vx) / G if vs else None
    sl = f"{-s:.2f}" if s else "  -"
    if ceil is None:
        v = "no Vswim in main text"
    elif z - 2*sem > ceil:
        v = f"IMPOSSIBLE: holo height above ceiling by {(z-ceil)/sem:.1f} SEM"
    elif z > ceil:
        v = "above ceiling, within 2 SEM"
    else:
        v = "consistent"
    cs = f"{ceil:.2f}" if ceil else "  -"
    print(f"{n:15s} {z:5.2f}+/-{sem:4.2f} {sl:>8s} {cs:>8s}  {v}")
