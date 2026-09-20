"""THE CONTROL.  SST_SG-G is (box mean) minus (global mean).  Since the box
term turns out to be near zero, the index's trend is essentially -1 x the
global trend -- which every box on Earth also gets.  So: how unusual is
-0.55 degC/century?  Slide a box of the same size over the whole ocean and
look at where the real subpolar gyre falls in that distribution.
"""
import numpy as np, hadisst
from ersst import trend

h = hadisst.HadISST(); WINTER=[11,12,1,2,3]
yg, g = h.season(h.box(-90,90,0,360), WINTER)
mg = (yg>=1871)&(yg<=2024)
gtr = trend(yg[mg], g[mg])[0]
print(f"global Nov-Mar SST trend 1871-2024: {gtr:+.3f} degC/century")

DLAT, DLON = 15, 35          # the SG box is 46-61N (15 deg) x 55-20W (35 deg)
res = []
for s in range(-60, 61, 5):
    for w in range(-180, 180, 10):
        v = h.box(s, s+DLAT, w, w+DLON)
        if np.ma.getmaskarray(v).any(): continue
        y, a = h.season(np.ma.filled(v, np.nan), WINTER)
        m = (y>=1871)&(y<=2024)
        if not np.isfinite(a[m]).all(): continue
        res.append((trend(y[m], a[m])[0] - gtr, s, w))
res.sort()
t = np.array([r[0] for r in res])
print(f"\n{len(res)} ocean boxes of 15 x 35 deg, index = boxmean - globalmean, 1871-2024")
print(f"  median {np.median(t):+.3f}   5th {np.percentile(t,5):+.3f}   95th {np.percentile(t,95):+.3f}  degC/century")
sg = -0.546
print(f"\n  real SST_SG-G = {sg:+.3f}")
print(f"  fraction of ocean boxes with an index trend at least as negative: "
      f"{(t <= sg).mean()*100:.1f}%  ({int((t<=sg).sum())} of {len(t)})")
print("\nfive most negative boxes:")
for x,s,w in res[:5]:
    print(f"  {x:+7.3f}  {s:3d}-{s+DLAT:3d}N  {w:4d}-{w+DLON:4d}E")
np.save('null_boxes.npy', np.array(res))
