import numpy as np, traj
from traj import atm as atm0
best=None
for f in (0.8,0.9,1.0,1.1,1.2,1.3):
  traj.atm=lambda z,f=f:(lambda p,r,T,c:(p*f,r*f,T,c))(*atm0(z))
  for fpa in (-7.8,-8.0,-8.2,-8.4,-8.6):
    s=traj.run(fpa0=fpa); v,_,z=s.y
    v62=np.interp(62e3,z[::-1],v[::-1]); v44=np.interp(44e3,z[::-1],v[::-1])
    e=abs(v62/11.2e3-1)+abs(v44/2.6e3-1)
    if best is None or e<best[0]: best=(e,f,fpa,v62,v44)
print(best)
