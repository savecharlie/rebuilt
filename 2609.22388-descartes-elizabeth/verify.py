"""Verification of arXiv:2609.22388 (Parrochia) — Descartes' Circle Theorem,
Princess Elizabeth, and Spinors.

Instrument: sympy exact rational/symbolic algebra.
LAW 1: the instrument is validated against answers known independently BEFORE
it is allowed to say anything new.
"""
import sympy as sp

d, e, f, x = sp.symbols('d e f x', positive=True)

# ---------------------------------------------------------------- GATE 0
# Known answer #1: the classic Apollonian quadruple (-1, 2, 2, 3).
k = [-1, 2, 2, 3]
assert sum(k)**2 == 2*sum(ki**2 for ki in k), "gate 0a"

# Known answer #2: three unit circles mutually tangent; inner Soddy curvature
# is 3 + 2*sqrt(3) (hand-derived: u^2 - 6u - 3 = 0).
u = sp.symbols('u')
sol = sp.solve(sp.Eq(3 + u**2, 6 + 6*u), u)
assert sp.simplify(max(sol) - (3 + 2*sp.sqrt(3))) == 0, "gate 0b"
print("GATE 0  known Apollonian answers reproduced          OK")

# ---------------------------------------------------------------- GATE 1
# Descartes eq. (4) as printed in the paper == Coxeter eq. (6) == the modern
# (sum k)^2 = 2 sum k^2, for the radii parametrisation k = 1/r.
eq4_L = d**2*e**2*f**2 + d**2*e**2*x**2 + d**2*f**2*x**2 + e**2*f**2*x**2
eq4_R = (2*d*e*f**2*x**2 + 2*d*e**2*f**2*x + 2*d*e**2*f*x**2
         + 2*d**2*e*f**2*x + 2*d**2*e*f*x**2 + 2*d**2*e**2*f*x)
eq4 = sp.expand(eq4_L - eq4_R)

ks = [1/d, 1/e, 1/f, 1/x]
modern = sp.expand(sp.together(sum(ks)**2 - 2*sum(ki**2 for ki in ks)) * (d*e*f*x)**2)
# modern numerator should be proportional to eq4
ratio = sp.simplify(sp.expand(modern) / eq4)
assert ratio.is_number, f"gate 1: not proportional, ratio={ratio}"
print(f"GATE 1  paper eq.(4) == modern Descartes, factor {ratio}   OK")

# ---------------------------------------------------------------- GATE 2
# A concrete numeric instance, built from COORDINATES, not from the formula.
# Three mutually tangent circles r=2,3,4; find the inner Soddy circle by
# solving the distance conditions directly, then confirm eq.(4) vanishes.
import itertools, math
rd, re_, rf = sp.Integer(2), sp.Integer(3), sp.Integer(4)
A = sp.Matrix([0, 0])
C = sp.Matrix([rd + rf, 0])
# B: |AB| = rd+re, |CB| = rf+re
bx = ((rd+re_)**2 + (rd+rf)**2 - (rf+re_)**2) / (2*(rd+rf))
by = sp.sqrt((rd+re_)**2 - bx**2)
B = sp.Matrix([bx, by])
hx, hy, hr = sp.symbols('hx hy hr', real=True)
sols = sp.solve([ (hx-A[0])**2 + (hy-A[1])**2 - (rd+hr)**2,
                  (hx-B[0])**2 + (hy-B[1])**2 - (re_+hr)**2,
                  (hx-C[0])**2 + (hy-C[1])**2 - (rf+hr)**2 ],
                [hx, hy, hr], dict=True)
inner = [s for s in sols if s[hr].is_real and s[hr] > 0]
assert len(inner) == 1, f"gate 2: expected one inner circle, got {sols}"
xr = sp.nsimplify(inner[0][hr])
chk = eq4.subs({d: rd, e: re_, f: rf, x: xr})
assert sp.simplify(chk) == 0, f"gate 2: eq(4) != 0 at coordinate solution ({chk})"
print(f"GATE 2  coordinate-built inner circle r={xr} satisfies eq.(4)  OK")
print()
print("="*66)
print("THE CLAIM: Parrochia's reconstruction of Descartes' missing step.")
print("="*66)

AC, AB, BC = d + f, d + e, e + f          # mutually tangent outer circles
HA, HB, HC = d + x, e + x, f + x          # H = inner circle centre, radius x

# Euclid II.13 feet of perpendiculars onto line AC
AK = sp.together((HA**2 + AC**2 - HC**2) / (2*AC))   # foot from H
AD = sp.together((AB**2 + AC**2 - BC**2) / (2*AC))   # foot from B

print("\nSTEP A  the two closed forms Descartes states in letter 2")
AK_paper = (d**2 + d*f + d*x - f*x)/(d + f)
AD_paper = (d**2 + d*f + d*e - f*e)/(d + f)
print("   AK  :", "OK" if sp.simplify(AK - AK_paper) == 0 else f"WRONG  {sp.simplify(AK-AK_paper)}")
print("   AD  :", "OK" if sp.simplify(AD - AD_paper) == 0 else f"WRONG  {sp.simplify(AD-AD_paper)}")

# heights above AC
HK2 = sp.expand(HA**2 - AK**2)            # = HC^2 - KC^2 too; check that
KC  = AC - AK
print("   HK^2 from A-side == HK^2 from C-side :",
      "OK" if sp.simplify(HK2 - (HC**2 - KC**2)) == 0 else "WRONG")
BD2 = sp.expand(AB**2 - AD**2)

KD = AD - AK

print("\nSTEP B  the Pythagorean closure  KD^2 + (BD-HK)^2 = HB^2")
Omega = sp.together(KD**2 + BD2 + HK2 - HB**2)      # = 2*BD*HK
print("   Omega (paper's simplified form)      :",
      "OK" if sp.simplify(Omega - 2*(d**2 + d*e + d*x - e*x - AK*AD)) == 0 else "WRONG")

Om_paper = (8*d*f*(d+f)*(e+x) - 8*(d**2+f**2)*e*x) / (d+f)**2
print("   Omega == paper's final fraction      :",
      "OK" if sp.simplify(Omega - Om_paper) == 0 else
      f"WRONG  residual {sp.simplify(sp.together(Omega-Om_paper))}")

N_paper = 4*d**2*f*e + 4*d**2*f*x + 4*d*f**2*e + 4*d*f**2*x - 4*d**2*e*x - 4*f**2*e*x
print("   N_Omega as printed                   :",
      "OK" if sp.simplify(sp.together(Omega) - 2*N_paper/(d+f)**2) == 0 else
      f"WRONG  {sp.factor(sp.simplify(sp.together(Omega)*(d+f)**2/2 - N_paper))}")

print("\nSTEP C  THE STEP THE PAPER WRITES WITH ELLIPSES:")
print("        square both sides, 4*BD^2*HK^2 == Omega^2, and reduce.")
residual = sp.expand(sp.together(4*BD2*HK2 - Omega**2) * (d+f)**4)
residual = sp.simplify(sp.expand(residual))
print("   residual factors as:")
print("      ", sp.factor(residual))
q = sp.simplify(sp.expand(residual) / eq4)
print("   residual / (paper's eq.4)  =", sp.factor(sp.simplify(q)))
