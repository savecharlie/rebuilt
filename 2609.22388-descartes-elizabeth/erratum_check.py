import sympy as sp
d,e,f,x = sp.symbols('d e f x', positive=True)
AC,AB,BC = d+f, d+e, e+f
HA,HB,HC = d+x, e+x, f+x
AK = (HA**2+AC**2-HC**2)/(2*AC); AD = (AB**2+AC**2-BC**2)/(2*AC)
HK2 = sp.expand(HA**2-AK**2); BD2 = sp.expand(AB**2-AD**2); KD = AD-AK
eq4 = sp.expand(d**2*e**2*f**2 + d**2*e**2*x**2 + d**2*f**2*x**2 + e**2*f**2*x**2
      - 2*d*e*f**2*x**2 - 2*d*e**2*f**2*x - 2*d*e**2*f*x**2
      - 2*d**2*e*f**2*x - 2*d**2*e*f*x**2 - 2*d**2*e**2*f*x)

Om_true  = sp.together(KD**2 + BD2 + HK2 - HB**2)
Om_paper = (8*d*f*(d+f)*(e+x) - 8*(d**2+f**2)*e*x)/(d+f)**2

print("Omega_paper / Omega_true              =", sp.simplify(Om_paper/Om_true))
print("Omega_true == 2*BD*HK (sign-checked)  =",
      sp.simplify(sp.expand(Om_true**2 - 4*BD2*HK2) / eq4))
print()
print("CLOSURE with the CORRECT Omega:")
r_true = sp.simplify(sp.expand((4*BD2*HK2 - Om_true**2)*(d+f)**4))
print("   (4 BD^2 HK^2 - Omega^2)*(d+f)^4  =", sp.factor(r_true))
print("   == -16 (d+f)^2 * [eq.4] ?          ", sp.simplify(r_true + 16*(d+f)**2*eq4) == 0)
print()
print("CLOSURE with the PAPER's PRINTED Omega:")
r_pap = sp.simplify(sp.expand((4*BD2*HK2 - Om_paper**2)*(d+f)**4))
print("   factors as                         ", sp.factor(r_pap))
print("   proportional to eq.4 ?             ", (sp.simplify(r_pap/eq4)).free_symbols == set() )
print()
# the HK^2 = f^2 - KC^2 typo
KC = AC-AK
print("paper step-4 'HK^2 = f^2 - KC^2' is off by:",
      sp.factor(sp.simplify(sp.expand(HK2 - (f**2 - KC**2)))))
