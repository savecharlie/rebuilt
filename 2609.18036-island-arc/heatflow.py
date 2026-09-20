"""Does 10 TW of earthquake-induced heat show up at the subduction zones, as the paper says it does?

Data:
  GHFDB Release 2024 v.2026.03 (GFZ Data Services, doi:10.5880/fidgeo.2024.014, CC-BY-4.0)
  GCMT catalogue 1976-2020 (globalcmt.org) -- used to LOCATE the subduction zones by the paper's
  own criterion: thrust events are those with Mrr > 0.
"""
import numpy as np
from parse import read_ndk

# ---------- heat flow ----------
rows = [l.rstrip("\n").split("\t") for l in
        open("ghfdb/IHFC_2024_GHFDB_v.2026.03.txt", encoding="utf-8", errors="replace")
        if not l.startswith("#")]
hdr = rows[5]                      # row 6 of the non-comment block is the column-name header
data = rows[6:]
ix = {n: i for i, n in enumerate(hdr)}
def col(r, n):
    v = r[ix[n]].strip() if ix[n] < len(r) else ""
    return v

q, lat, lon, env, elev = [], [], [], [], []
for r in data:
    try:
        qq = float(col(r, "q")); la = float(col(r, "lat_NS")); lo = float(col(r, "long_EW"))
    except ValueError:
        continue
    if not np.isfinite(qq) or abs(la) > 90 or abs(lo) > 180:
        continue
    try:    ev = float(col(r, "elevation"))
    except ValueError: ev = np.nan
    q.append(qq); lat.append(la); lon.append(lo); env.append(col(r, "environment")); elev.append(ev)

q = np.array(q); lat = np.array(lat); lon = np.array(lon)
env = np.array(env, dtype=object); elev = np.array(elev)
keep = (q > -50) & (q < 1000)
q, lat, lon, env, elev = q[keep], lat[keep], lon[keep], env[keep], elev[keep]
print(f"heat-flow points used: {len(q)}   median q = {np.median(q):.0f} mW/m^2")

# ---------- thrust earthquakes (the paper's own locator for subduction) ----------
e   = read_ndk("jan76_dec20.ndk")
M0  = e["M0"] * 1e-7
Mw  = (2/3) * (np.log10(M0) - 9.1)
thr = (e["Mrr"] > 0) & (Mw >= 6.0)
tlat, tlon = e["lat"][thr], e["lon"][thr]
print(f"thrust events Mw>=6 used as the subduction locator: {thr.sum()}")

# ---------- great-circle distance to the nearest thrust event ----------
R = 6371.0
def to_xyz(la, lo):
    la, lo = np.radians(la), np.radians(lo)
    return np.stack([np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo), np.sin(la)], -1)
A, B = to_xyz(lat, lon), to_xyz(tlat, tlon)
dmin = np.empty(len(A))
for i in range(0, len(A), 2000):                      # chunk: 91k x 5k would be 3.6 GB
    c = np.clip(A[i:i+2000] @ B.T, -1, 1)
    dmin[i:i+2000] = R * np.arccos(c).min(1)

# ---------- is the near-trench heat flow anomalous? ----------
offsh = np.array([("offshore" in str(s)) or ("marine" in str(s).lower()) for s in env])
print(f"\ndistance to nearest thrust Mw>=6      all points          offshore only")
print(f"   (km)                        N   median  mean      N   median  mean")
for lo_, hi in [(0,100),(100,200),(200,300),(300,500),(500,1000),(1000,3000),(3000,1e9)]:
    m = (dmin >= lo_) & (dmin < hi)
    mo = m & offsh
    lab = f"{lo_:>5.0f}-{hi:<6.0f}" if hi < 1e8 else f"{lo_:>5.0f}+     "
    s1 = f"{m.sum():6d} {np.median(q[m]):7.0f} {q[m].mean():6.0f}" if m.sum() else " "*21
    s2 = f"{mo.sum():6d} {np.median(q[mo]):7.0f} {q[mo].mean():6.0f}" if mo.sum() else " "*21
    print(f"  {lab}   {s1}   {s2}")
np.save("dmin.npy", dmin); np.save("q.npy", q); np.save("offsh.npy", offsh)
