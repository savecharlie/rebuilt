"""What the authors' own mapping says, computed exactly instead of argued.
Every number here comes from the gated instrument in ising.py."""
import numpy as np, math, json
from ising import Cluster, geometry

EPS_B, D_EPS_B, D_EPS, EPS_A = -3.5, 2.5, 4.0, 3.0
h_of_c = lambda c: (math.log(c) - EPS_B - D_EPS_B) / 2.0
c_of_h = lambda h: math.exp(2 * h + EPS_B + D_EPS_B)
CSTAR = c_of_h(0.0)
cl = Cluster(5); pin = cl.pin_centre_up()
OUT = {}

print("="*78)
print("0.  THE FIELD AXIS.  h = (ln c - eps_b - d_eps_b)/2 = (ln c + 1)/2")
print("="*78)
print(f"  zero effective field at c* = e^-1 = {CSTAR:.6f}")
print(f"  their headline point  c = 0.32, J_b = 0.55  ->  h = {h_of_c(0.32):+.6f}")
print("  the whole c-range of their Fig. 4 in field units:")
for c in (0.10, 0.15, 0.20, 0.25, 0.32, 0.368, 0.45, 0.60):
    print(f"      c = {c:.3f}  ->  h = {h_of_c(c):+.4f}")
print(f"  ... so their ENTIRE explored concentration axis spans "
      f"dh = {h_of_c(0.60)-h_of_c(0.10):.3f} in field.")
OUT["c_star"] = CSTAR; OUT["h_at_headline"] = h_of_c(0.32)

print()
print("="*78)
print("1.  IS THE UNACTIVATED, LIGAND-FREE STATE S0 THE EQUILIBRIUM STATE? (P1)")
print("="*78)
print("   J_b      c      h       <bound>/N  (no seed)     <bound>/N  (one correct ligand)")
grid = {}
for J in (0.30, 0.40, 0.45, 0.50, 0.55, 0.60, 0.70):
    for c in (0.20, 0.32, 0.368, 0.45):
        h = h_of_c(c)
        a, b = cl.mean_bound(J, h), pin.mean_bound(J, h)
        grid[(J, c)] = (a, b)
        if c in (0.32, 0.368):
            print(f"  {J:5.2f}  {c:5.3f}  {h:+.4f}      {a:.6f}                {b:.6f}")
hp = h_of_c(0.32)
OUT["headline_noseed"] = cl.mean_bound(0.55, hp)
OUT["headline_seeded"] = pin.mean_bound(0.55, hp)
print(f"\n  AT THEIR HEADLINE POINT (J_b=0.55, c=0.32):")
print(f"    no correct ligand : <n_bound> = {25*OUT['headline_noseed']:.4f} of 25  "
      f"({OUT['headline_noseed']:.2%})")
print(f"    one correct ligand: <n_bound> = {25*OUT['headline_seeded']:.4f} of 25  "
      f"({OUT['headline_seeded']:.2%})")
print(f"    the seed IS one of those, so it recruits "
      f"{25*OUT['headline_seeded']-1:.3f} beyond itself "
      f"(and raises the total by {25*(OUT['headline_seeded']-OUT['headline_noseed']):.3f})")
print(f"    (their simulations report the cluster reaching Na = 10.)")

print()
print("="*78)
print("2.  WHERE IS THE EQUILIBRIUM FLIP FIELD?  (P3)")
print("="*78)
print("  h_half: the field at which the seed-free cluster is half bound.  If their")
print("  c_max sat at equilibrium this is where it would be.")
def h_half(J, obj, lo=-6.0, hi=6.0):
    for _ in range(200):
        m = (lo + hi) / 2
        if obj.mean_bound(J, m) < 0.5: lo = m
        else: hi = m
    return (lo + hi) / 2
print("    J_b     h_half      c_half      (c* = 0.3679)    h_half with a seed")
rows=[]
for J in (0.0, 0.20, 0.30, 0.40, 0.4407, 0.50, 0.55, 0.60, 0.70, 0.85, 1.00):
    hh, hs = h_half(J, cl), h_half(J, pin)
    rows.append((J, hh, c_of_h(hh), hs))
    print(f"   {J:6.3f}  {hh:+.5f}   {c_of_h(hh):.5f}                       {hs:+.5f}")
OUT["h_half"] = rows

print()
print("="*78)
print("3.  THE EXACT NUCLEATION BARRIER vs THEIR HAND-ARGUED 4*J_b  (P4)")
print("="*78)
print("  At a = 0 there is NO barrier: F(n) climbs monotonically, so a scalar dF+ is")
print("  meaningless (argmax returns n=25 always, which is just 40J - 50h).  The right")
print("  objects are the STEPS.  Both of their hand-argued barriers, checked exactly:")
print("    J_b   F(1)-F(0)   8J-2h-ln25    first step beside seed   -ln[4e^-(4J-2h)+20e^-(8J-2h)]")
for J in (0.40, 0.55, 0.70, 0.85):
    F = cl.Fn(J, hp); Fs = pin.Fn(J, hp)
    iso = 8*J - 2*hp - math.log(25)
    grow = -math.log(4*math.exp(-(4*J-2*hp)) + 20*math.exp(-(8*J-2*hp)))
    print(f"   {J:5.2f}   {F[1]-F[0]:9.6f}  {iso:11.6f}    {Fs[2]-Fs[1]:14.6f}   {grow:20.6f}")
print("\n  and the correct ligand's exact advantage in reaching a bound decamer -> 8*J_b:")
print("    J_b   F10 no seed   F10 seeded   advantage   /J_b    rate factor")
adv_tab=[]
for J in (0.40, 0.50, 0.55, 0.60, 0.70, 0.85, 1.00):
    Fu = cl.Fn(J, hp); Fs = pin.Fn(J, hp)
    u = Fu[10]-np.nanmin(Fu); v = Fs[10]-np.nanmin(Fs); a = u-v
    adv_tab.append((J,u,v,a))
    print(f"   {J:5.2f}   {u:10.4f}   {v:10.4f}   {a:9.4f}  {a/J:6.3f}   x{math.exp(a):8.1f}")
OUT["seed_advantage"] = adv_tab
print()
print("="*78)
print("4.  WHAT ACTUALLY HOLDS THE AMPLIFIED STATE UP: THE ACTIVATION RATCHET")
print("="*78)
print("  Table I: r_off(a) = r_off(0) exp(-d_eps * a), d_eps = 4.  Detailed balance at")
print("  fixed a then gives K_d(a) = exp(eps_b + d_eps_b - d_eps*a), so")
print("       h(a) = h(0) + d_eps*a/2 .")
for a, name in ((0.0,"a=0   (unactivated)"), (0.5,"a=1/2 (one step)"), (1.0,"a=1   (signalling)")):
    hh = hp + D_EPS*a/2
    print(f"    {name}:  r_off = {0.37*math.exp(-D_EPS*a):.4f}/t0   h = {hh:+.4f}   "
          f"c-equivalent = {c_of_h(hh):8.4f}   <bound>/N = {cl.mean_bound(0.55,hh):.6f}")
span = h_of_c(0.50)-h_of_c(0.15)          # the c axis of their Fig. 4
print(f"\n  h(a) - h(0) = d_eps*a/2:  one half-step = {D_EPS*0.5/2:.1f}, "
      f"full activation = {D_EPS*1.0/2:.1f}.")
print(f"  Their Fig. 4 concentration axis (c = 0.15 to 0.50) is worth {span:.3f} in field,")
print(f"  so full activation outruns the whole axis by {D_EPS*1.0/2/span:.2f}x and even one")
print(f"  half-step by {D_EPS*0.5/2/span:.2f}x.")
OUT["ratchet_field_per_full_activation"] = D_EPS/2

print()
print("="*78)
print("5.  WHY 5x5.  THE CORRELATION LENGTH AT THEIR OPERATING COUPLING")
print("="*78)
for J in (0.40, 0.4407, 0.50, 0.55, 0.60, 0.70, 0.85):
    Ks = 0.5*math.asinh(1/math.sinh(2*J)) if J > 0 else float('inf')
    xi = 1/(2*(J-Ks)) if J > 0.4407 else float('inf')
    print(f"    J_b = {J:.4f}   xi = {xi:8.4f} lattice sites = {xi*10:8.2f} nm"
          + ("   <-- their operating point" if abs(J-0.55)<1e-9 else ""))
Ks=0.5*math.asinh(1/math.sinh(2*0.55)); xi55=1/(2*(0.55-Ks))
print(f"\n  cluster is 5 sites = 50 nm across; xi(J_b=0.55) = {xi55:.3f} sites = {xi55*10:.1f} nm")
print(f"  -> the cluster is {5/xi55:.2f} correlation lengths wide.")
print(f"  (Their own ref. [5]: TCR preclusters are 30-300 nm.)")
OUT["xi_at_0p55"] = xi55
json.dump(OUT, open("RESULTS.json","w"), indent=1, default=float)
print("\nwrote RESULTS.json")
