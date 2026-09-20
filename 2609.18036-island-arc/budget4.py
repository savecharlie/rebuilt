"""Second instrument. The first (budget3.py) failed its own control: integrating the heat-flow
excess around thrust and around normal seismicity gave 6.0 and 6.2 TW, indistinguishable, when
the paper puts 10 TW at one and -2 TW at the other. That integral was measuring the
continent/ocean contrast, not a plate-boundary anomaly.

So: oceanic points only, and normalise each one against the heat flow its own seafloor AGE
predicts, with age taken from the water depth the database already records (GDH1). Then the
only thing left is the residual, and thrust vs normal is an internally controlled comparison -
same database, same normalisation, same everything but the sign of Mrr.

Known bias, stated because it runs the wrong way for the paper's critics: a trench is deep by
flexure, not by age, so depth->age reads too OLD there, predicts too LITTLE heat, and inflates
the residual near subduction. Any deficit found is therefore conservative.
"""
import numpy as np
from scipy.spatial import cKDTree
from parse import read_ndk

Rk = 6371.0
def xyz(la, lo):
    la, lo = np.radians(la), np.radians(lo)
    return np.stack([np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo), np.sin(la)], -1)*Rk

rows = [l.rstrip("\n").split("\t") for l in
        open("ghfdb/IHFC_2024_GHFDB_v.2026.03.txt", encoding="utf-8", errors="replace")
        if not l.startswith("#")]
hdr, data = rows[5], rows[6:]
ix = {n: i for i, n in enumerate(hdr)}
Q=[];LAT=[];LON=[];EL=[]
for r in data:
    try:
        qq=float(r[ix["q"]]); la=float(r[ix["lat_NS"]]); lo=float(r[ix["long_EW"]]); el=float(r[ix["elevation"]])
    except (ValueError, IndexError): continue
    if -50 < qq < 1000 and abs(la)<=90 and abs(lo)<=180 and np.isfinite(el):
        Q.append(qq);LAT.append(la);LON.append(lo);EL.append(el)
Q,LAT,LON,EL = map(np.array,(Q,LAT,LON,EL))

# --- oceanic only, and age from depth (GDH1: d = 2600 + 365 sqrt(t), t < 20 Ma;
#     d = 5651 - 2473 exp(-0.0278 t) beyond) ---
oc = EL < -2600
d  = -EL[oc]
t_young = ((d - 2600.0)/365.0)**2
with np.errstate(invalid="ignore"):
    arg = (5651.0 - d)/2473.0
t_old = np.where(arg > 0, -np.log(np.clip(arg, 1e-9, None))/0.0278, np.nan)
age = np.where(t_young < 20.0, t_young, t_old)
ok  = np.isfinite(age) & (age > 0.5) & (age < 180)
# GDH1 heat flow
q_pred = np.where(age < 20.0, 510.0/np.sqrt(np.clip(age,0.5,None)),
                  48.0 + 96.0*np.exp(-0.0278*age))
qo, lato, lono = Q[oc][ok], LAT[oc][ok], LON[oc][ok]
res = qo - q_pred[ok]
print(f"oceanic points with a usable depth->age: {ok.sum()}")
print(f"  median observed {np.median(qo):.0f}, median predicted {np.median(q_pred[ok]):.0f}, "
      f"median residual {np.median(res):+.0f} mW/m^2")

e = read_ndk("jan76_dec20.ndk"); M0=e["M0"]*1e-7; Mw=(2/3)*(np.log10(M0)-9.1)
P = xyz(lato, lono)
out = {}
for name, mask in (("thrust  (Mrr>0)", (e["Mrr"]>0)&(Mw>=6.0)),
                   ("normal  (Mrr<0)", (e["Mrr"]<0)&(Mw>=6.0))):
    T = cKDTree(xyz(e["lat"][mask], e["lon"][mask]))
    c,_ = T.query(P, workers=-1)
    out[name] = 2*Rk*np.arcsin(np.clip(c/(2*Rk),-1,1))

print("\n  distance      thrust: N   median residual      normal: N   median residual")
for a,b in [(0,50),(50,100),(100,200),(200,300),(300,500),(500,1000),(1000,99999)]:
    line=f"  {a:5d}-{b if b<9e4 else 0:<5d}"
    for name in out:
        m=(out[name]>=a)&(out[name]<b)
        line += f"   {m.sum():7d} {np.median(res[m]):+8.0f}" if m.sum()>=30 else f"   {m.sum():7d}    (few)"
    print(line)
np.save("res.npy",res); np.save("d_thr.npy",out["thrust  (Mrr>0)"]); np.save("d_nor.npy",out["normal  (Mrr<0)"])
