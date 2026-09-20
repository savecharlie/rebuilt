"""Where, in depth, is the thrust seismic moment? And how much gravitational power is available
from slab descent over that depth range?"""
import numpy as np
from parse import read_ndk

e = read_ndk("jan76_dec20.ndk"); M0 = e["M0"]*1e-7
Mw = (2/3)*(np.log10(M0)-9.1)
thr = e["Mrr"] > 0
d, m = e["dep"][thr], M0[thr]
tot = m.sum()
print("thrust (Mrr>0) seismic moment by centroid depth, GCMT 1976-2020")
print("   depth (km)      events      % of thrust moment     cumulative")
cum = 0
for a,b in [(0,30),(30,50),(50,70),(70,100),(100,200),(200,400),(400,700),(700,900)]:
    s = (d>=a)&(d<b); f = 100*m[s].sum()/tot; cum += f
    print(f"   {a:4d}-{b:<5d}   {s.sum():7d}        {f:8.2f}            {cum:7.2f}")
print(f"\n  moment-weighted mean depth: {(d*m).sum()/tot:.1f} km")
for z in (70, 100, 200):
    print(f"  fraction of thrust moment shallower than {z:3d} km: {100*m[d<z].sum()/tot:.1f}%")

# --- gravitational power of slab descent ---
print("\n--- gravitational power available from slab descent ---")
A_dot = 3.0e6          # m^2/yr of plate subducted (= seafloor creation rate, Parsons 1981 /
                       #   Rowley 2002 give 3.0-3.4 km^2/yr)
h_slab = 100e3         # m, thermal lithosphere
drho   = 50.0          # kg/m^3, ~500 K colder at alpha = 3e-5, rho = 3300
g      = 9.8
YR     = 365.25*86400
Vdot   = A_dot*h_slab                      # m^3/yr
for H_km in (70, 100, 660, 2890):
    P = drho*g*(H_km*1e3)*Vdot/YR
    print(f"  descent through the top {H_km:5d} km: {P/1e12:6.2f} TW")
print("  (the 2890 km figure is the whole-mantle number; Morgan et al. 2016, as the paper quotes"
      " them, give 11-14 TW for the descending limb)")

# --- what island-arc volcanism actually needs ---
print("\n--- heat required by island-arc volcanism ---")
rho_m, cp, L, dT = 2800.0, 1200.0, 4.0e5, 1200.0
h_per_m3 = rho_m*(cp*dT + L)
for V_km3 in (2.0, 3.0, 6.0, 8.0):
    P = V_km3*1e9*h_per_m3/YR
    print(f"  {V_km3:4.1f} km^3/yr of arc magma: {P/1e12:5.2f} TW"
          f"   ({h_per_m3/1e9:.2f} GJ per m^3: {cp*dT/1e3:.0f} kJ/kg sensible + {L/1e3:.0f} latent)")
