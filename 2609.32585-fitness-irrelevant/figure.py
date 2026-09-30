import numpy as np, sys, json
sys.path.insert(0,'.')
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from switching import Switching
from ceiling import ceiling_composition

P  = lambda x: 0.03*np.cos(4*np.pi*x)+0.03*x
r0 = lambda x: 0.1+1.9/(1+np.exp(-(x-0.5)/0.02))
g  = lambda U: 1.0-U; gp = lambda U: -1.0
BG="#f4efe6"; INK="#22201d"; ACC="#b4472e"; BLU="#2f5d7c"; GRY="#9a938a"

# ---- (a) the tilt --------------------------------------------------------
s = Switching(1.0,300,lambda x:0.02+0*x,P=P); rv=r0(s.x)
G = np.log(1/0.01)
pc,sstar = ceiling_composition(s.psi,rv,s.h,G)
t,Y = s.run(0.01*s.psi,r0,g,T=40,n_out=1600,gp=gp,rtol=1e-11)
U=s.mass(Y); p=Y/U[:,None]; d=s.h*np.abs(p-s.psi).sum(axis=1); k=d.argmax()

# ---- (b) the climb -------------------------------------------------------
climb=[]
for D0 in (0.4,0.2,0.1,0.05,0.03,0.02,0.012,0.008):
    sd=Switching(1.0,300,lambda x:D0+0*x,P=P); lam,_=sd.spectral_gap()
    rvd=r0(sd.x); rb=sd.h*np.sum(rvd*sd.psi)
    pcd,_=ceiling_composition(sd.psi,rvd,sd.h,G); cl=sd.h*np.abs(pcd-sd.psi).sum()
    tt,YY=sd.run(0.01*sd.psi,r0,g,T=60,n_out=2000,gp=gp,rtol=1e-8)
    UU=sd.mass(YY); pp=YY/UU[:,None]; dd=sd.h*np.abs(pp-sd.psi).sum(axis=1)
    climb.append((rb/lam, dd.max()/cl)); print('  climb D=%g eps=%.2f frac=%.4f'%(D0,rb/lam,dd.max()/cl),flush=True)
climb=np.array(climb)

# ---- (c) the sign --------------------------------------------------------
sign=[]
rbar_psi = s.h*np.sum(rv*s.psi)
for U0 in (0.003,0.01,0.03,0.1,0.3,0.6,0.85,1.2,1.7,3.0,8.0,30.0):
    G2=np.log(1.0/U0); pc2,ss2=ceiling_composition(s.psi,rv,s.h,G2)
    tt,YY=s.run(U0*s.psi,r0,g,T=60,n_out=2000,gp=gp,rtol=1e-8)
    UU=s.mass(YY); pp=YY/UU[:,None]
    shift=s.h*(pp@rv)-rbar_psi
    j=np.abs(shift).argmax()
    sign.append((G2, shift[j], s.h*np.sum(rv*pc2)-rbar_psi)); print('  sign U0=%g G=%.3f meas=%+.4f pred=%+.4f'%(U0,G2,shift[j],s.h*np.sum(rv*pc2)-rbar_psi),flush=True)
sign=np.array(sign)
np.savez("figure_data.npz",climb=climb,sign=sign,pc=pc,peak=p[k],psi=s.psi,x=s.x,
         sstar=sstar,tpeak=t[k])

fig,ax=plt.subplots(1,3,figsize=(14.5,4.6),facecolor=BG)
for a in ax:
    a.set_facecolor(BG)
    for sp in ("top","right"): a.spines[sp].set_visible(False)
    for sp in ("left","bottom"): a.spines[sp].set_color(GRY)
    a.tick_params(colors=INK,labelsize=9)

A=ax[0]
A.fill_between(s.x,0,s.psi,color=BLU,alpha=.13,lw=0)
A.plot(s.x,s.psi,color=BLU,lw=1.8)
A.plot(s.x,pc,color=ACC,lw=3.2,alpha=.45)
A.plot(s.x,p[k],color=INK,lw=1.4,ls=(0,(4,2.2)))
A.set_xlabel("phenotype  $x$",color=INK,fontsize=10)
A.set_ylabel("composition  $p = u/U$",color=INK,fontsize=10)
A.text(.06,.93,"where selection gets to",transform=A.transAxes,fontsize=12.5,
       color=INK,style="italic",va="top")
A.text(.13,3.05,"$\\psi$",color=BLU,fontsize=13)
A.text(.60,5.4,"$\\psi\\,e^{r s^*}/M(s^*)$",color=ACC,fontsize=11.5)
A.text(.60,4.55,"measured peak, $t=%.2f$"%t[k],color=INK,fontsize=9)
A.set_ylim(0,None)

B=ax[1]
B.axhline(1,color=GRY,lw=.9,ls=":")
B.plot(climb[:,0],climb[:,1],"-",color=INK,lw=1.4,zorder=2)
B.plot(climb[:,0],climb[:,1],"o",ms=5,color=BG,mec=INK,mew=1.4,zorder=3)
i=np.argmin(np.abs(climb[:,0]-10.76))
B.plot(climb[i,0],climb[i,1],"o",ms=8,color=ACC,zorder=4)
B.annotate("the paper's own example\nsits at %.0f%% of the ceiling"%(100*climb[i,1]),
           (climb[i,0],climb[i,1]),textcoords="offset points",xytext=(-6,-52),
           fontsize=9,color=ACC,ha="center")
B.set_xscale("log"); B.set_ylim(0,1.12)
B.set_xlabel("$\\varepsilon=\\bar r_\\psi/\\lambda_1$   (selection speed / mixing speed)",
             color=INK,fontsize=10)
B.set_ylabel("peak imprint  /  ceiling",color=INK,fontsize=10)
B.text(.05,.93,"it climbs, and never crosses",transform=B.transAxes,fontsize=12.5,
       color=INK,style="italic",va="top")

C=ax[2]
C.axhline(0,color=GRY,lw=.9); C.axvline(0,color=GRY,lw=.9)
C.plot(sign[:,0],sign[:,2],"-",color=ACC,lw=3.2,alpha=.45)
C.plot(sign[:,0],sign[:,1],"o-",ms=4.5,color=INK,lw=1.3,mfc=BG,mew=1.2)
C.set_xlabel("$\\ln(U^*/U_0)$   —   e-foldings of growth",color=INK,fontsize=10)
C.set_ylabel("shift in mean fitness at the peak",color=INK,fontsize=10)
C.text(.05,.93,"a shrinking population\nselects backwards",transform=C.transAxes,
       fontsize=12.5,color=INK,style="italic",va="top")
C.text(-4.6,-.30,"shrinks  →  enriched for\nthe SLOWEST phenotypes",fontsize=8.6,color=INK)
C.text(1.2,.75,"grows  →  enriched\nfor the fastest",fontsize=8.6,color=INK)
C.text(3.4,-.28,"ceiling",color=ACC,fontsize=9)

fig.suptitle("Selection's whole budget is the growth"
             "        $K(s^*)=\\ln(U^*/U_0)$,   $K$ = cumulant generating function of fitness under $\\psi$",
             color=INK,fontsize=12.5,y=1.005)
fig.text(.995,-.02,"Iris · fire 310 · rebuilt from Fassoni arXiv:2609.32585",
         ha="right",color=GRY,fontsize=8)
fig.tight_layout()
fig.savefig("selection_budget.png",dpi=155,facecolor=BG,bbox_inches="tight")
print("s* =",sstar,"peak t =",t[k])
print("climb:",np.round(climb,4).tolist())
print("sign:",np.round(sign,4).tolist())
