import numpy as np
from scipy.spatial import cKDTree
from parse import read_ndk

Rk, Rm = 6371.0, 6371e3
def xyz(la, lo):
    la, lo = np.radians(la), np.radians(lo)
    return np.stack([np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo), np.sin(la)], -1)*Rk
def chord(r): return 2*Rk*np.sin(r/(2*Rk))

# --- heat-flow points ---
rows = [l.rstrip("\n").split("\t") for l in
        open("ghfdb/IHFC_2024_GHFDB_v.2026.03.txt", encoding="utf-8", errors="replace")
        if not l.startswith("#")]
hdr, data = rows[5], rows[6:]
ix = {n: i for i, n in enumerate(hdr)}
Q, LAT, LON = [], [], []
for r in data:
    try:
        qq=float(r[ix["q"]]); la=float(r[ix["lat_NS"]]); lo=float(r[ix["long_EW"]])
    except (ValueError, IndexError): continue
    if -50 < qq < 1000 and abs(la) <= 90 and abs(lo) <= 180:
        Q.append(qq); LAT.append(la); LON.append(lo)
Q, LAT, LON = map(np.array, (Q, LAT, LON))
P = xyz(LAT, LON)

# --- equal-area grid for band areas ---
nlon, nlat = 2880, 1440
lonc = np.linspace(-180,180,nlon,endpoint=False)+180/nlon
latc = np.degrees(np.arcsin(np.linspace(-1,1,nlat,endpoint=False)+1.0/nlat))
cell = 4*np.pi*Rm**2/(nlon*nlat)
LO, LA = np.meshgrid(lonc, latc)
G = xyz(LA.ravel(), LO.ravel())

e = read_ndk("jan76_dec20.ndk"); M0 = e["M0"]*1e-7
Mw = (2/3)*(np.log10(M0)-9.1)
sets = {"THRUST (subduction), paper: 10 TW": ((e["Mrr"]>0)&(Mw>=6.0), 10.0),
        "NORMAL (spreading),  paper:  2 TW": ((e["Mrr"]<0)&(Mw>=6.0),  2.0)}

EDGES = [0,50,100,150,200,300,400,500,700,1000,1500]
for label,(mask,claim) in sets.items():
    T = cKDTree(xyz(e["lat"][mask], e["lon"][mask]))
    dq,_ = T.query(P, workers=-1);  dq = 2*Rk*np.arcsin(np.clip(dq/(2*Rk),-1,1))
    dg,_ = T.query(G, workers=-1);  dg = 2*Rk*np.arcsin(np.clip(dg/(2*Rk),-1,1))
    base = np.median(Q[dq > 2000])
    print(f"\n=== {label} ===")
    print(f"  baseline: median q at >2000 km from any such event = {base:.0f} mW/m^2 (N={(dq>2000).sum()})")
    print("   annulus (km)    N    median q   excess    area 10^12 m^2    excess x area")
    tot = 0.0
    for a,b in zip(EDGES[:-1], EDGES[1:]):
        m  = (dq>=a)&(dq<b)
        ar = ((dg>=a)&(dg<b)).sum()*cell
        if m.sum() < 30:
            print(f"   {a:5d}-{b:<5d} {m.sum():6d}        (too few)"); continue
        exc = np.median(Q[m]) - base
        P_TW = exc*1e-3*ar/1e12
        tot += max(P_TW, 0.0)
        print(f"   {a:5d}-{b:<5d} {m.sum():6d} {np.median(Q[m]):9.0f} {exc:+8.0f}   {ar/1e12:10.1f}"
              f"      {P_TW:+8.2f} TW")
    print(f"   integrated positive excess out to {EDGES[-1]} km: {tot:.2f} TW   "
          f"(the paper needs {claim:.0f} TW)")
