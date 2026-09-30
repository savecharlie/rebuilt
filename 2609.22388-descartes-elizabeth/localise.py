import sympy as sp
d,e,f,x = sp.symbols('d e f x', positive=True)

AK = (d**2 + d*f + d*x - f*x)/(d+f)
AD = (d**2 + d*f + d*e - f*e)/(d+f)

# Paper, step 1: Omega = 2*( d^2+de+dx-ex - AK*AD ), everything over (d+f)^2
block   = d**2 + d*e + d*x - e*x
N_true  = sp.expand(block*(d+f)**2 - (d**2+d*f+d*x-f*x)*(d**2+d*f+d*e-f*e))

print("TRUE numerator N_Omega  (block*(d+f)^2  -  num(AK)*num(AD)):")
print("   ", sp.factor(N_true), "\n    expanded:", N_true)

N_paper = 4*d**2*f*e + 4*d**2*f*x + 4*d*f**2*e + 4*d*f**2*x - 4*d**2*e*x - 4*f**2*e*x
print("\nPAPER's printed N_Omega:")
print("   ", N_paper)
print("\n   N_paper / N_true =", sp.simplify(N_paper/N_true))

# Now check the paper's own two printed expansions, step 2 and step 3.
s2_true  = sp.expand(block*(d**2 + 2*d*f + f**2))
s2_paper = (d**4 + 2*d**3*f + d**2*f**2 + d**3*e + 2*d**2*f*e + d*f**2*e
            + d**3*x + 2*d**2*f*x + d*f**2*x - d**2*e*x - 2*d*f*e*x - f**2*e*x)
print("\nSTEP 2 (printed 12 monomials):",
      "OK" if sp.expand(s2_true-s2_paper)==0 else f"WRONG by {sp.expand(s2_true-s2_paper)}")

s3_true  = sp.expand(-(d**2+d*f+d*x-f*x)*(d**2+d*f+d*e-f*e))
s3_paper = (-d**4 - 2*d**3*f - d**2*f**2 - d**3*e
            - d**2*f*e + d**2*f*e + d*f**2*e - d**3*x - d**2*f*x + d**2*f*x
            + d*f**2*x - d**2*e*x + 2*d*f*e*x - f**2*e*x)
print("STEP 3 (printed 12 monomials):",
      "OK" if sp.expand(s3_true-s3_paper)==0 else f"WRONG by {sp.expand(s3_true-s3_paper)}")
print("   sum of the two PRINTED expansions =", sp.expand(s2_paper+s3_paper))
print("   ... which is N_true? ", sp.expand(s2_paper+s3_paper-N_true)==0)
