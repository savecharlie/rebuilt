"""Rebuild k_i and alpha_i per station from public tables (Silber & Bowman 2025
Tables S2/S3) + reconstructed trajectory, then test the §4.2 'Fourth' inference."""
import numpy as np
from scipy.stats import spearmanr
from traj import run, atm
def load(fn, cols):
    d={}
    for ln in open(fn):
        p=ln.split(); d[p[0]]=[float(p[c]) for c in cols]
    return d
s2=load('s2b.txt',[5,6]); s3=load('s3b.txt',[1,5])
st=sorted(set(s2)&set(s3)); print(len(st),'stations')
tau=np.array([s2[k][0] for k in st]); dtau=np.array([s2[k][1] for k in st])
zs=np.array([s3[k][0] for k in st]); R=np.array([s3[k][1] for k in st])
s=run(); v,_,z=s.y
M=np.array([np.interp(zz*1e3, z[::-1], (v/np.array([atm(q)[3] for q in z]))[::-1]) for zz in zs])
R0=23.99*tau-0.73          # the paper's Eq.12 — reproduces its inversion to 2.3% MAPE
d=0.81; alpha=R0/(M*d); k=alpha*d/np.sqrt(0.7*1.0*0.515)
print("median alpha %.3f (paper 0.294)  median k %.3f (paper 0.396)"%(np.median(alpha),np.median(k)))
print("spearman k vs z: %.2f (paper -0.84)"%spearmanr(k,zs)[0])
print("tau range %.3f-%.3f  z range %.1f-%.1f  R range %.0f-%.0f km"%(tau.min(),tau.max(),zs.min(),zs.max(),R.min(),R.max()))
# how much of the k trend is tau vs M?
print("spearman tau vs z %.2f, M vs z %.2f"%(spearmanr(tau,zs)[0],spearmanr(M,zs)[0]))
p=np.polyfit(np.log(tau),np.log(M),1)[0]
print("d lnM / d ln tau along trajectory: %.2f   (inversion's d lnR0/d ln tau = 1.35 per paper)"%p)
for pp in (1.35,p):
    kk=tau**pp/M; print(" if R0 ∝ tau^%.2f: spearman(k,z)=%.2f  spread(IQR/med)=%.2f"%(pp,spearmanr(kk,zs)[0],np.subtract(*np.percentile(kk,[75,25]))/np.median(kk)))
np.savez('stations.npz',st=st,tau=tau,zs=zs,R=R,M=M,R0=R0,k=k,alpha=alpha)
