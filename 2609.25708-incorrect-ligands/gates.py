"""Gates.  The instrument is not allowed to say anything new until every one is green.
Law 1 (known answer first), law 6 (two independent routes), law 1c (scan the knob)."""
import numpy as np, itertools, math, sys
from ising import (Cluster, geometry, transfer_logZ, site_fields,
                   onsager_logZ_per_site, strip_logZ_per_site, enumerate_histogram)

FAILS = []
def gate(name, ok, detail):
    print(("  PASS  " if ok else "  FAIL  ") + name + " :: " + detail)
    if not ok: FAILS.append(name)

def brute_logZ(L, J, h, h_seed=0.0):
    """Third, dumbest route: build every state explicitly and sum e^-E from the
    Hamiltonian as written, with no B/P decomposition anywhere."""
    bonds, miss = geometry(L); n = L*L; centre = (L//2)*L + (L//2)
    tot = -np.inf
    for bits in itertools.product([-1,1], repeat=n):
        s = np.array(bits)
        E = -J*sum(s[i]*s[j] for i,j in bonds)
        for i in range(n):
            E -= (h - J*miss[i] + (h_seed if i==centre else 0.0)) * s[i]
        tot = np.logaddexp(tot, -E)
    return float(tot)

print("== gate 1: enumeration bookkeeping ==")
H5, _ = enumerate_histogram(5)
gate("all 2^25 states accounted for", H5.sum() == 1<<25, f"{H5.sum()} == {1<<25}")
cl = Cluster(5)

print("== gate 2: J = 0 must be independent sites (exact analytic) ==")
for h in (-1.0, -0.0697, 0.0, 0.5):
    got = cl.mean_bound(0.0, h)
    want = (1 + math.tanh(h)) / 2
    gate(f"J=0, h={h:+.4f}", abs(got-want) < 1e-12, f"{got:.15f} vs tanh form {want:.15f}")

print("== gate 3: two independent routes -- histogram vs transfer matrix (law 6) ==")
for (J, h, hs) in [(0.0,0.3,0.0),(0.55,-0.0697,0.0),(0.55,-0.0697,5.0),
                   (0.30,0.20,2.5),(0.70,-0.40,-3.0),(1.2,0.05,0.0)]:
    a = cl.logZ(J,h,hs); b = transfer_logZ(5,J,site_fields(5,J,h,hs))
    gate(f"logZ J={J} h={h} hseed={hs}", abs(a-b) < 1e-9, f"hist {a:.12f} vs TM {b:.12f}  d={a-b:+.2e}")

print("== gate 4: third route -- explicit 3x3 Hamiltonian, no decomposition ==")
c3 = Cluster(3)
for (J,h,hs) in [(0.55,-0.0697,0.0),(0.4,0.3,2.5),(0.8,-0.5,-1.0)]:
    a = c3.logZ(J,h,hs); b = brute_logZ(3,J,h,hs); c = transfer_logZ(3,J,site_fields(3,J,h,hs))
    gate(f"3x3 three routes J={J} h={h} hseed={hs}",
         abs(a-b)<1e-10 and abs(a-c)<1e-10, f"hist {a:.12f} brute {b:.12f} TM {c:.12f}")

print("== gate 5: saturation limits ==")
gate("h -> +inf all bound", abs(cl.mean_bound(0.55, 40.0) - 1.0) < 1e-12, f"{cl.mean_bound(0.55,40.0):.15f}")
gate("h -> -inf all free",  abs(cl.mean_bound(0.55,-40.0) - 0.0) < 1e-12, f"{cl.mean_bound(0.55,-40.0):.15e}")
gate("J=0,h=0 half bound",  abs(cl.mean_bound(0.0, 0.0) - 0.5) < 1e-14, f"{cl.mean_bound(0.0,0.0):.15f}")

print("== gate 6: pinned centre really is a pinned centre ==")
p = cl.pin_centre_up()
gate("pin == h_seed -> +inf", abs(p.mean_bound(0.55,-0.0697) - cl.mean_bound(0.55,-0.0697,60.0)) < 1e-12,
     f"{p.mean_bound(0.55,-0.0697):.15f} vs {cl.mean_bound(0.55,-0.0697,60.0):.15f}")

print("== gate 7: EXTERNAL known answer -- Onsager 1944 ==")
Kc = math.log(1+math.sqrt(2))/2
gate("critical coupling sinh(2Kc)=1", abs(math.sinh(2*Kc)-1) < 1e-14, f"Kc={Kc:.12f}")
# the famous exact internal energy at T_c:  u = -sqrt(2) J
d = 1e-5
u = -(onsager_logZ_per_site(Kc+d) - onsager_logZ_per_site(Kc-d))/(2*d)
gate("Onsager u(K_c) = -sqrt(2)", abs(u + math.sqrt(2)) < 3e-3, f"{u:.6f} vs {-math.sqrt(2):.6f}")
print("     cylinder -> bulk.  Law 4: the reference is a suspect too -- two candidate closed")
print("     forms for xi, and the DATA picks which.  Law 1c: the regime number is W/xi; scan it.")
print("     SEAM (law 10): over W in [4,12] the prefactor exponent p and xi are strongly")
print("     correlated in the fit, so p is NOT determined here and no value is claimed for it.")
Ws=np.array([4,6,8,10,12]); rows=[]
for K in (0.50, 0.55, 0.60, 0.70, 0.85, 1.00):
    ons=onsager_logZ_per_site(K, 4000)
    d=np.array([strip_logZ_per_site(int(W),None,K)[0] for W in Ws])-ons
    A=np.c_[np.ones(len(Ws)), -np.log(Ws), -Ws]           # d = A W^-p e^-W/xi, p FREE
    coef,*_=np.linalg.lstsq(A, np.log(d), rcond=None)
    xi=1/coef[2]; res=np.max(np.abs(np.exp(A@coef)/d-1))
    Ks=0.5*math.asinh(1/math.sinh(2*K))                    # dual: sinh2K sinh2K* = 1
    c2, c4 = 1/(2*(K-Ks)), 1/(4*(K-Ks))
    rows.append(dict(K=K, xi=xi, c2=c2, c4=c4, res=res, dl=d[-1], reg=Ws[-1]/c2, p=coef[1]))
    print(f"       K={K:.2f}  W/xi={Ws[-1]/c2:5.2f}  xi={xi:.4f} vs exact {c2:.4f} "
          f"(ratio {xi/c2:.4f})   cand4 {c4:.4f} (ratio {xi/c4:.3f})   resid {res:.3%}  p~{coef[1]:.2f}")
rows.sort(key=lambda r: r["reg"])
gate("3-parameter form describes every K to better than 0.5%",
     all(r["res"] < 0.005 for r in rows), "max %.3f%%" % (100*max(r["res"] for r in rows)))
gate("measured xi -> exact 1/[2(K-K*)] as the regime number W/xi grows",
     all(rows[i]["xi"]/rows[i]["c2"] > rows[i+1]["xi"]/rows[i+1]["c2"] for i in range(len(rows)-1))
     and abs(rows[-1]["xi"]/rows[-1]["c2"]-1) < 0.005,
     " -> ".join(f'{r["xi"]/r["c2"]:.4f}' for r in rows) + f'  (deepest: {abs(rows[-1]["xi"]/rows[-1]["c2"]-1):.2%})')
gate("rival form 1/[4(K-K*)] rejected by ~2x everywhere, drifting further away",
     all(abs(r["xi"]/r["c4"]-1) > 0.4 for r in rows),
     " ".join(f'{r["xi"]/r["c4"]:.3f}' for r in rows))
gate("cylinder free energy -> Onsager from above, monotonically in W/xi",
     all(r["dl"] > 0 for r in rows) and all(rows[i]["dl"] > rows[i+1]["dl"] for i in range(len(rows)-1)),
     " > ".join(f'{r["dl"]:+.1e}' for r in rows))

print("\n== gate 8: law 1c -- scan the knob the gates fixed ==")
worst = 0.0
for J in np.arange(0.0, 1.51, 0.1):
    for h in (-2.0,-0.5,-0.0697,0.0,0.5,2.0):
        d_ = abs(cl.logZ(J,h) - transfer_logZ(5,J,site_fields(5,J,h)))
        worst = max(worst, d_)
gate("hist==TM over the whole (J,h) grid actually used", worst < 1e-8, f"worst |dlogZ| = {worst:.2e}")

print()
if FAILS:
    print("RED:", FAILS); sys.exit(1)
print("ALL GATES GREEN -- the instrument may now say something new.")
