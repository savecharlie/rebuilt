"""Section 3 arithmetic of arXiv:2609.21714, recomputed. Iris, fire 290."""
import numpy as np
G=6.67430e-11
print("P_imp(Naiad, >10 km) = 3.3e-7 x 3.5e-5 =", f"{3.3e-7*3.5e-5:.3g} /yr  (text: 1.2e-11)")
D=2*33e3; v=20.7e3
for lab,Q in (("asteroid-like 6e8 erg/g",6e4),("KBO-like 3e7 erg/g",3e3)):
    print(f"D_imp {lab}: {(2*Q/v**2)**(1/3)*D/1e3:.2f} km  (text: 4 / 1.5)")
rhoP=6835100e9/G/(4/3*np.pi*24764e3**3); delta=48227/24764*(800/rhoP)**(1/3)
f=1-(1.5/delta)**3
print(f"rho_Nep={rhoP:.0f}; delta={delta:.3f} (1.54); Q*_TD/Q*_D={f:.3f} -> factor {1/f:.1f} (~10);"
      f" KBO impactor {(2*3e3*f/v**2)**(1/3)*D/1e3:.2f} km (~0.6)")
print(f"q_imp = 0.8(q+1), q=-1.7+-0.2 -> {0.8*(-1.7+1):.2f} +- {0.8*0.2:.2f}  (text: -0.56 +- 0.16)")
P2=1.155e-11*(0.2)**-2.1
for q in (-0.4,-0.56,-0.72):
    for dk in (1.5,0.6):
        print(f"  q_imp={q}: mean time to a >{dk} km impactor = {1/(P2*(dk/2)**q)/1e9:.2f} Gyr")
