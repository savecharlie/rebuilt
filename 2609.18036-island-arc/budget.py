"""Global moment/energy budget from GCMT 1976-2020, against the paper's three numbers:
   seismic power ~10 GW, thrust-induced Eg release ~10 TW, normal-faulting Eg gain ~2 TW."""
import numpy as np
from parse import read_ndk

DC = 1e-7                      # dyne-cm -> N m
YR = 365.25 * 86400.0
e   = read_ndk("jan76_dec20.ndk")
M0  = e["M0"] * DC
Mrr = e["Mrr"] * DC
yrs = 45.0                     # 1976/01/01 .. 2020/12/31

thr = Mrr > 0                  # the paper's own criterion: positive Mrr = thrust
nor = Mrr < 0

print(f"GCMT 1976-2020: {len(e)} events, {yrs:.0f} yr")
print(f"  thrust (Mrr>0) {thr.sum():6d}   normal (Mrr<0) {nor.sum():6d}")
print()
tot_rate = M0.sum() / yrs
print(f"total moment rate        {tot_rate:.3e} N m / yr")
print(f"  thrust                 {M0[thr].sum()/yrs:.3e}   ({100*M0[thr].sum()/M0.sum():.0f}%)")
print(f"  normal                 {M0[nor].sum()/yrs:.3e}   ({100*M0[nor].sum()/M0.sum():.0f}%)")
print()

# --- seismic energy: Kanamori's E_s = M0 / (2e4) ---
Es_rate = M0.sum() / 2e4 / yrs           # J/yr
print(f"seismic energy (Kanamori E_s = M0/2e4):  {Es_rate/YR/1e9:.2f} GW"
      f"    [paper says ~10 GW]")

# --- Eg, two independent calibrations of |dEg|/M0 ---
# (a) the paper's own Sumatra figure: 3e21 J for the 2004/12/26 event
sel = e[e["date"] == "2004/12/26"]; sum_M0 = sel["M0"].max() * DC
r_sumatra = 3e21 / sum_M0
# (b) Dahlen 1977 as the paper quotes it: |dEg| ~ 1000 x E_s, and E_s = M0/2e4
r_dahlen  = 1000.0 / 2e4
print(f"\n|dEg|/M0 calibrations:  Sumatra {r_sumatra:.4f}   Dahlen-1000x {r_dahlen:.4f}")

for name, r in (("Sumatra-calibrated", r_sumatra), ("Dahlen 1000x", r_dahlen)):
    for lab, m in (("thrust  (Eg release)", M0[thr].sum()), ("normal  (Eg gain)", M0[nor].sum())):
        P = r * m / yrs / YR
        print(f"  {name:18s}  {lab:22s} {P/1e12:7.3f} TW")
print("\n  [paper says ~10 TW release / ~2 TW gain]")
