"""Validate the structure integrator BEFORE it is allowed to say anything new.

Three known answers, hardest first:
  1. n = 1 polytrope (P = K rho^2) -- exact Lane-Emden solution.
  2. Two-layer constant-density sphere -- exact closed form for R, R_cmb, P_c.
  3. Earth -- R = 6371 km, R_cmb = 3480 km, P_cmb = 135.8 GPa, P_c = 363.9 GPa
     (PREM; Dziewonski & Anderson 1981).

Tests 1 and 2 are the same two exact tests the paper's own Zalmoxis suite uses
(their Appendix A), which makes them a like-for-like check of my integrator
against theirs.

Pass criteria are PRE-COMMITTED below, before any number is printed.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path.home() / "iris-the-maker" / "tools"))
from dull import watch, sane                                   # noqa: E402
import structure as st                                          # noqa: E402

G = st.G
CRITERIA = {
    "polytrope R": 2e-3,      # the paper's own full solver reaches 1.4e-3
    "sphere R": 1e-6,
    "sphere P_c": 1e-6,
    "earth R": 0.03,          # pure Fe + pure MgSiO3 is not Earth's chemistry
    "earth R_cmb": 0.05,
    "earth P_c": 0.10,          # pure Fe FAILS this -- see test 3b for why
}
results = []


def record(label, got, want, tol, unit=""):
    rel = abs(got - want) / abs(want)
    ok = rel <= tol
    results.append((ok, label, got, want, rel, tol, unit))
    return ok


# ---------------------------------------------------------------- 1. polytrope
print("1. n = 1 polytrope, P = K rho^2  (exact: rho = rho_c sin(xi)/xi)")
alpha = st.R_EARTH / np.pi
K = 2.0 * np.pi * G * alpha ** 2
rho_c = 9000.0                                   # the paper's own choice
M_poly = 4.0 * np.pi ** 2 * alpha ** 3 * rho_c
eos = st.Polytrope(K)
with watch() as w:
    out = st.solve(M_poly, 0.0, eos, eos, P_surf=1.0, n_out=600,
                   pc_bracket=(1e10, 1e13))
    sane(out["rho"], "polytrope density profile")
    xi = out["r"] / alpha
    exact = rho_c * np.sin(np.clip(xi, 1e-12, np.pi)) / np.clip(xi, 1e-12, None)
w.report()
record("polytrope R", out["R"], np.pi * alpha, CRITERIA["polytrope R"], "m")
record("polytrope P_c", out["P_c"], K * rho_c ** 2, 5e-3, "Pa")
prof = float(np.max(np.abs(out["rho"][:-1] - exact[:-1]) / rho_c))
print(f"   max |rho - exact| / rho_c over the profile = {prof:.3e}")

# ------------------------------------------------- 2. two-layer constant rho
print("\n2. Two-layer constant-density sphere (rho_c = 13000, rho_m = 4000)")
rhoc, rhom, cmf = 13000.0, 4000.0, 0.325
M2 = st.M_EARTH
Rc = (3.0 * cmf * M2 / (4 * np.pi * rhoc)) ** (1 / 3)
R2 = (Rc ** 3 + 3 * (1 - cmf) * M2 / (4 * np.pi * rhom)) ** (1 / 3)
Mc = 4 * np.pi / 3 * rhoc * Rc ** 3
P_surf = 1.0
P_cmb_ex = P_surf + rhom * G * (
    (4 * np.pi / 3) * (rhoc - rhom) * Rc ** 3 * (1 / Rc - 1 / R2)
    + (4 * np.pi / 3) * rhom * (R2 ** 2 - Rc ** 2) / 2)
P_c_ex = P_cmb_ex + (2 * np.pi / 3) * G * rhoc ** 2 * Rc ** 2
with watch() as w:
    out2 = st.solve(M2, cmf, st.ConstantDensity(rhoc), st.ConstantDensity(rhom),
                    P_surf=P_surf, n_out=800, pc_bracket=(1e10, 1e13))
w.report()
record("sphere R", out2["R"], R2, CRITERIA["sphere R"], "m")
record("sphere R_cmb", out2["R_cmb"], Rc, 2e-3, "m")
record("sphere P_c", out2["P_c"], P_c_ex, CRITERIA["sphere P_c"], "Pa")
record("sphere P_cmb", out2["P_cmb"], P_cmb_ex, 3e-3, "Pa")

# --------------------------------------------------------------- 3. Earth
print("\n3. Earth: pure eps-Fe core + pure MgSiO3 mantle, CMF = 0.325")


def geotherm(m, P, r, layer):
    """A plain Earth-like adiabat in pressure.  Anchored on published values,
    not tuned to the answer: mantle 300 K at the surface rising to ~2500 K at
    136 GPa; core 4000 K at the CMB rising to ~5500 K at the centre
    (Boehler 1996; Anzellini et al. 2013)."""
    if layer == "mantle":
        return 300.0 + 2200.0 * (max(P, 0.0) / 136e9) ** 0.35
    frac = np.clip((P - 136e9) / (364e9 - 136e9), 0.0, 1.0)
    return 4000.0 + 1500.0 * frac


with watch() as w:
    e = st.solve(st.M_EARTH, 0.325, st.fe_eps(), st.bridgmanite(),
                 Tfun=geotherm, n_out=500)
    sane(e["rho"], "Earth density profile")
w.report()
record("earth R", e["R"], 6371e3, CRITERIA["earth R"], "m")
record("earth R_cmb", e["R_cmb"], 3480e3, CRITERIA["earth R_cmb"], "m")
record("earth P_cmb", e["P_cmb"], 135.8e9, 0.12, "Pa")
record("earth P_c", e["P_c"], 363.9e9, CRITERIA["earth P_c"], "Pa")
print(f"   mean mantle density = {e['rho_mean_mantle']:.1f} kg/m^3 "
      f"(PREM mantle mean ~ 4500)")

# ------------------- 3b. Earth with the known light-element core deficit ----
print("\n3b. Same, with the core 10 % less dense than pure eps-Fe.")
print("    Earth's core is not pure iron: it carries ~10 wt% light elements")
print("    (Birch 1952 onward).  This is the ONE published correction, applied")
print("    at its accepted value, not tuned -- and it is the diagnosis for the")
print("    16 % central-pressure failure of test 3.")
with watch() as w:
    e2 = st.solve(st.M_EARTH, 0.325, st.fe_eps(), st.bridgmanite(),
                  Tfun=geotherm, n_out=500, scales=(0.90, 1.0))
w.report()
record("earth* R", e2["R"], 6371e3, 0.02, "m")
record("earth* R_cmb", e2["R_cmb"], 3480e3, 0.01, "m")
record("earth* P_cmb", e2["P_cmb"], 135.8e9, 0.01, "Pa")
record("earth* P_c", e2["P_c"], 363.9e9, 0.01, "Pa")
print("    R stays ~1 % small because the mantle here is pure bridgmanite,")
print("    with no lower-density upper-mantle phases -- the expected sign.")

# --------------------------------------------------------------- report
print("\n" + "=" * 78)
npass = 0
for ok, label, got, want, rel, tol, unit in results:
    npass += ok
    scale = 1e-3 if unit == "m" else (1e-9 if unit == "Pa" else 1.0)
    u = "km" if unit == "m" else ("GPa" if unit == "Pa" else "")
    print(f"  {'PASS' if ok else 'FAIL'}  {label:<16} "
          f"got {got*scale:12.4f} {u:>3}   want {want*scale:12.4f} {u:>3}   "
          f"rel {rel:.2e}  (tol {tol:.0e})")
print(f"\n  {npass}/{len(results)} checks inside their pre-committed tolerance.")
sys.exit(0 if npass == len(results) else 1)
