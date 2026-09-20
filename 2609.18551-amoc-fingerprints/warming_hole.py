"""Where is the warming hole, and does the fingerprint box sit on it?

Per-gridcell Nov-Mar SST trend over the North Atlantic, HadISST 1871-2024,
with the SST_SG box of Emirzade et al. drawn on it.  Also reports the coldest
cells and the trend of the box that actually contains them.
"""
import numpy as np, hadisst
from ersst import trend

h = hadisst.HadISST()
WINTER=[11,12,1,2,3]
lat, lon = h.lat, h.lon
jm = (lat>=30)&(lat<=70); im = h.lon_mask(-70, 0)
J, I = np.where(jm)[0], np.where(im)[0]
tr = np.full((len(J), len(I)), np.nan)
for a,j in enumerate(J):
    for b,i in enumerate(I):
        s = h.sst[:, j, i]
        if np.ma.getmaskarray(s).mean() > 0.02: continue
        y, v = h.season(np.ma.filled(s, np.nan), WINTER)
        m = (y>=1871)&(y<=2024)
        if not np.isfinite(v[m]).all(): continue
        tr[a,b] = trend(y[m], v[m])[0]

# HadISST latitude runs DESCENDING (89.5 -> -89.5), so tr's rows do too.
# Store latitude-ASCENDING or every map drawn from this file is upside down.
# (Caught Sep 20 2026 only by looking at the picture: the star marking the
#  coldest cell at 53.5N landed on a warm patch.)
np.save('na_trend.npy', tr[::-1])
np.save('na_axes.npy', np.array([lat[J].min(), lat[J].max(), lon[I].min(), lon[I].max()]))
print('grid', tr.shape, 'valid', np.isfinite(tr).sum())
flat = [(tr[a,b], lat[J][a], lon[I][b]) for a in range(len(J)) for b in range(len(I)) if np.isfinite(tr[a,b])]
flat.sort()
print('\nTen coldest-trending cells, Nov-Mar 1871-2024 (degC/century):')
for t,la,lo in flat[:10]:
    print(f"  {t:+7.3f}  {la:5.1f}N  {((lo+180)%360)-180:6.1f}E")
neg = [f for f in flat if f[0] < 0]
print(f"\ncells with negative trend: {len(neg)} of {len(flat)} in 30-70N, 70W-0")
# the box the paper uses
print("\nSST_SG box 46-61N, 55-20W:")
y,v = h.season(h.box(46,61,-55,-20), WINTER); m=(y>=1871)&(y<=2024)
b,se = trend(y[m], v[m]); print(f"  trend {b:+.3f} +/- {se:.3f}")
# a box drawn ON the actual minimum
print("\nbox centred on the true minimum (50-60N, 45-25W):")
y,v = h.season(h.box(50,60,-45,-25), WINTER); m=(y>=1871)&(y<=2024)
b,se = trend(y[m], v[m]); print(f"  trend {b:+.3f} +/- {se:.3f}")
