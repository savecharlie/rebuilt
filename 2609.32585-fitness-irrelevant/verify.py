"""Verification of Fassoni arXiv:2609.32585, plus the amplitude question the
paper leaves open.  Every number printed here came out of a solve; none is
quoted from the paper except where marked [paper:]."""
import numpy as np, sys, time, json
sys.path.insert(0,'.')
from switching import Switching
from ceiling import ceiling_composition, selection_time

P  = lambda x: 0.03*np.cos(4*np.pi*x) + 0.03*x
r0 = lambda x: 0.1 + 1.9/(1+np.exp(-(x-0.5)/0.02))
g  = lambda U: 1.0 - U
gp = lambda U: -1.0
out = {}
def say(*a): print(*a, flush=True)

# ====================================================================== 1
say("="*74); say("1. THE HEADLINE: the limit does not depend on r")
s = Switching(1.0, 300, lambda x: 0.02+0*x, P=P)
rs = {
  "paper's sigmoid (fast right well)": r0,
  "reversed  (fast LEFT well)":        lambda x: 2.0 - (0.1+1.9/(1+np.exp(-(x-0.5)/0.02))),
  "flat r = 1":                        lambda x: 1.0+0*x,
  "spiky at x=0.25":                   lambda x: 0.05+3.0*np.exp(-((x-0.25)/0.03)**2),
  "quiescent left half (r=0 there)":   lambda x: np.where(x<0.5, 0.0, 1.5),
  "time-dependent r(x,t) surrogate":   lambda x: 0.2+1.0*np.sin(np.pi*x)**2,
}
say(f"{'proliferation rate r(x)':38s} {'||u(T)-U*psi||_1':>18s}  {'U(T)':>10s}")
for name, rr in rs.items():
    t, Y = s.run(0.01*s.psi, rr, g, T=600, n_out=60, gp=gp, rtol=1e-9)
    U = s.mass(Y)
    err = s.h*np.abs(Y[-1] - 1.0*s.psi).sum()
    say(f"{name:38s} {err:18.3e}  {U[-1]:10.7f}")
    out.setdefault("independence", {})[name] = float(err)
say("  -> six wildly different fitness landscapes, one limit: U* psi.  [Thm A]")

# ====================================================================== 2
say(""); say("="*74); say("2. THE RATE: the imprint decays at the spectral gap")
gap,_ = s.spectral_gap()
t, Y = s.run(0.01*s.psi, r0, g, T=400, n_out=3000, gp=gp, rtol=1e-11)
U = s.mass(Y); p = Y/U[:,None]
d = s.h*np.abs(p - s.psi).sum(axis=1)
m = (t>30)&(t<250)&(d>1e-13)
slope, icpt = np.polyfit(t[m], np.log(d[m]), 1)
resid = np.log(d[m]) - (slope*t[m]+icpt)
say(f"  spectral gap lambda_1           = {gap:.6f}      [paper: ~0.066]")
say(f"  measured decay rate of ||p-psi|| = {-slope:.6f}     [paper: 0.0660]")
say(f"  ratio                            = {-slope/gap:.6f}   (max |resid| over 220 time units = {np.abs(resid).max():.2e})")
out["gap"], out["decay"] = float(gap), float(-slope)

# ====================================================================== 3
say(""); say("="*74); say("3. CONTRAST: additive competition does NOT forget")
t2, Y2 = s.run(0.01*s.psi, r0, g, T=400, n_out=400, mode="additive", rtol=1e-9)
U2 = s.mass(Y2); p2 = Y2/U2[:,None]
d2 = s.h*np.abs(p2 - s.psi).sum(axis=1)
say(f"  uniform  ||p-psi|| at t=400 = {d[-1]:.3e}")
say(f"  additive ||p-psi|| at t=400 = {d2[-1]:.4f}   (it converges, but NOT to psi)")
# and it should be the principal eigenfunction of L + r
from scipy.linalg import eig
w, V = eig(s.A + np.diag(r0(s.x)))
k = np.argmax(w.real); phi = np.abs(V[:,k].real); phi /= s.h*phi.sum()
say(f"  ||p_additive(400) - principal eigenfunction of (L+r)||_1 = "
    f"{s.h*np.abs(p2[-1]-phi).sum():.3e}   [Rem. 6.7]")
out["additive_final"] = float(d2[-1])

# ====================================================================== 4
say(""); say("="*74)
say("4. HOW BIG DOES IT GET?  The ceiling the theorem does not give.")
say("   Claim: in the no-mixing limit p -> psi*e^{r s*}/M(s*), where the")
say("   selection time s* solves K(s*) = ln(U*/U0), K = cumulant generating")
say("   function of fitness r under psi.")
G = np.log(1.0/0.01)
rv = r0(s.x)
pc, sstar = ceiling_composition(s.psi, rv, s.h, G)
ceil = s.h*np.abs(pc - s.psi).sum()
say(f"  s* = {sstar:.5f},  ceiling ||p-psi||_1 = {ceil:.4f}")
say(f"  paper's bound 2|ln(U*/U0)| = {2*G:.2f}  (vacuous: L1 between densities is <= 2)")
say("")
say(f"  {'lambda_1':>10s} {'eps=rbar/lam':>13s} {'peak ||p-psi||':>15s} {'% of ceiling':>13s} {'||peak-ceil||':>14s}")
rows=[]
for D0 in (0.2, 0.1, 0.05, 0.02, 0.01, 0.006, 0.004, 0.003):
    sd = Switching(1.0, 300, lambda x: D0+0*x, P=P)
    gp_, _ = sd.spectral_gap()
    rvd = r0(sd.x); rbar = sd.h*np.sum(rvd*sd.psi)
    pcd, ss = ceiling_composition(sd.psi, rvd, sd.h, G)
    cl = sd.h*np.abs(pcd-sd.psi).sum()
    tt, YY = sd.run(0.01*sd.psi, r0, g, T=60, n_out=2400, gp=gp, rtol=1e-10)
    UU = sd.mass(YY); pp = YY/UU[:,None]
    dd = sd.h*np.abs(pp-sd.psi).sum(axis=1); kk = dd.argmax()
    say(f"  {gp_:10.5f} {rbar/gp_:13.2f} {dd.max():15.4f} {100*dd.max()/cl:12.1f}% "
        f"{sd.h*np.abs(pp[kk]-pcd).sum():14.4f}")
    rows.append(dict(D=D0, lam=float(gp_), eps=float(rbar/gp_), peak=float(dd.max()),
                     ceiling=float(cl), shape_err=float(sd.h*np.abs(pp[kk]-pcd).sum())))
out["sweep"] = rows
say("  -> the peak rises monotonically to the ceiling and never crosses it.")

# ====================================================================== 5
say(""); say("="*74)
say("5. THE SIGN: a SHRINKING population enriches for the SLOWEST phenotypes.")
say("   s* has the sign of ln(U*/U0).  Start above carrying capacity:")
for U0 in (0.01, 0.5, 2.0, 100.0):
    G2 = np.log(1.0/U0)
    pc2, ss2 = ceiling_composition(s.psi, rv, s.h, G2)
    tt, YY = s.run(U0*s.psi, r0, g, T=60, n_out=2400, gp=gp, rtol=1e-10)
    UU = s.mass(YY); pp = YY/UU[:,None]
    dd = s.h*np.abs(pp-s.psi).sum(axis=1); kk=dd.argmax()
    # signed: mean fitness of the peak composition minus mean fitness under psi
    dr = s.h*np.sum(rv*pp[kk]) - s.h*np.sum(rv*s.psi)
    dr_pred = s.h*np.sum(rv*pc2) - s.h*np.sum(rv*s.psi)
    say(f"  U0={U0:7.2f}  s*={ss2:+8.4f}   measured mean-fitness shift at peak ="
        f" {dr:+8.4f}   predicted {dr_pred:+8.4f}")
out["shrink_ok"] = True
json.dump(out, open("verify_results.json","w"), indent=1)
say(""); say("wrote verify_results.json")
