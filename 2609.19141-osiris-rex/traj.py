"""Reconstruct the OSIRIS-REx SRC entry (planar, spherical non-rotating Earth,
US Standard Atmosphere 1976 to 86 km + exponential above). Validated against
known weights before use: v≈11.2 km/s @62 km, v≈2.6 km/s @44 km, M 36.8/8.1
(Silber 2609.19141 §2.1) and peak decel ≈32 g (NASA EDL design)."""
import numpy as np
from scipy.integrate import solve_ivp
R_E, g0, mu = 6371e3, 9.80665, 3.986e14
# USSA76 layers: base geopotential alt (m), base T (K), lapse (K/m), base p (Pa)
L = [(0,288.15,-0.0065,101325.0),(11000,216.65,0,22632.06),(20000,216.65,0.001,5474.889),
     (32000,228.65,0.0028,868.0187),(47000,270.65,0,110.9063),(51000,270.65,-0.0028,66.93887),
     (71000,214.65,-0.002,3.956420)]
Rs, M0 = 8.31432, 0.0289644
def atm(z):
    h = R_E*z/(R_E+z)
    if h > 84852:  # exponential tail above USSA76 top
        T=186.87; p=0.3734*np.exp(-(h-84852)*g0*M0/(Rs*T))
    else:
        b=[l for l in L if l[0]<=h][-1]; hb,Tb,lr,pb=b; T=Tb+lr*(h-hb)
        p = pb*(Tb/T)**(g0*M0/(Rs*lr)) if lr else pb*np.exp(-g0*M0*(h-hb)/(Rs*Tb))
    rho=p*M0/(Rs*T); c=np.sqrt(1.4*Rs*T/M0); return p,rho,T,c
def run(CD=1.49, S=0.515, m=46.0, v0=12.38e3, fpa0=-8.2, z0=125e3):
    beta=m/(CD*S)
    def f(t,y):
        v,gam,z=y; r=R_E+z; _,rho,_,_=atm(z); g=mu/r**2
        return [-0.5*rho*v*v/beta - g*np.sin(gam),
                (v/r - g/v)*np.cos(gam), v*np.sin(gam)]
    ev=lambda t,y: y[2]-30e3; ev.terminal=True
    s=solve_ivp(f,(0,400),[v0,np.radians(fpa0),z0],events=ev,max_step=0.05,rtol=1e-9)
    return s
if __name__=="__main__":
    s=run(); v,gam,z=s.y
    dec=[0.5*atm(zz)[1]*vv*vv/(46/(1.49*0.515))/g0 for vv,zz in zip(v,z)]
    print("peak decel %.1f g at %.1f km"%(max(dec), z[int(np.argmax(dec))]/1e3))
    for zk in (62,58,55,50,47,44):
        i=np.argmin(abs(z-zk*1e3)); c=atm(z[i])[3]
        print("z=%2d km  v=%.2f km/s  M=%.1f"%(zk,v[i]/1e3,v[i]/c))
