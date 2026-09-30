import numpy as np, sys
sys.path.insert(0,'.')
from switching import Switching
from ceiling import ceiling_composition
from scipy.integrate import solve_ivp

P  = lambda x: 0.03*np.cos(4*np.pi*x)+0.03*x
r0 = lambda x: 0.1+1.9/(1+np.exp(-(x-0.5)/0.02))
g  = lambda U: 1.0-U
s  = Switching(1.0,300,lambda x:0.02+0*x,P=P)
gap,mu = s.spectral_gap()

print("A. IS THE TAIL A CLEAN EXPONENTIAL AT lambda_1?")
t,Y = s.run(0.01*s.psi, r0, g, T=400, n_out=4000, gp=lambda U:-1.0, rtol=1e-12)
U=s.mass(Y); p=Y/U[:,None]; d=s.h*np.abs(p-s.psi).sum(axis=1)
for lo,hi in [(30,250),(60,250),(100,250),(150,300),(200,350)]:
    m=(t>lo)&(t<hi)&(d>1e-13)
    sl,ic=np.polyfit(t[m],np.log(d[m]),1); res=np.log(d[m])-(sl*t[m]+ic)
    print(f"  window [{lo},{hi}]  rate={-sl:.6f}  rate/lambda_1={-sl/gap:.5f}  max|resid|={np.abs(res).max():.2e}")
print(f"  lambda_1={gap:.6f}   lambda_2={-np.sort(mu)[-3]:.6f}  (gap between them sets how fast the fit cleans up)")

print("\nB. THE CEILING, WITH SWITCHING TURNED COMPLETELY OFF (A=0).")
print("   Then the ODE is exactly dp/ds=(r-rbar)p, so p_inf must be psi*e^{r s*}/M(s*).")
rv=r0(s.x); h=s.h
for U0 in (0.01, 0.2, 2.0, 50.0):
    G=np.log(1.0/U0)
    pc,ss = ceiling_composition(s.psi,rv,h,G)
    def rhs(t,u): return g(h*u.sum())*rv*u          # NO switching at all
    sol=solve_ivp(rhs,(0,4000),U0*s.psi,method="BDF",rtol=1e-12,atol=1e-16,
                  jac=lambda t,u: np.diag(g(h*u.sum())*rv) - h*np.outer(rv*u,np.ones(len(u))))
    uf=sol.y[:,-1]; pf=uf/(h*uf.sum())
    print(f"  U0={U0:6.2f}  s*={ss:+9.4f}   ||p_inf(numeric) - psi e^{{r s*}}/M||_1 = "
          f"{h*np.abs(pf-pc).sum():.3e}    U_inf={h*uf.sum():.10f}")
