import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
d=np.load("figure_data.npz")
climb,sign,pc,peak,psi,x,sstar,tpk = (d["climb"],d["sign"],d["pc"],d["peak"],
                                      d["psi"],d["x"],float(d["sstar"]),float(d["tpeak"]))
BG="#f4efe6"; INK="#22201d"; ACC="#b4472e"; BLU="#2f5d7c"; GRY="#9a938a"
fig,ax=plt.subplots(1,3,figsize=(14.8,4.8),facecolor=BG)
for a in ax:
    a.set_facecolor(BG)
    for sp in ("top","right"): a.spines[sp].set_visible(False)
    for sp in ("left","bottom"): a.spines[sp].set_color(GRY)
    a.tick_params(colors=INK,labelsize=9)

A=ax[0]
A.fill_between(x,0,psi,color=BLU,alpha=.13,lw=0)
A.plot(x,psi,color=BLU,lw=1.9)
A.plot(x,pc,color=ACC,lw=4.0,alpha=.42,solid_capstyle="round")
A.plot(x,peak,color=INK,lw=1.3,ls=(0,(5,2.6)))
A.set_xlabel("phenotype  $x$",color=INK,fontsize=10)
A.set_ylabel("composition  $p=u/U$",color=INK,fontsize=10)
A.set_ylim(0,6.6); A.set_xlim(0,1)
A.text(.03,6.25,"where selection gets to",fontsize=12.5,color=INK,style="italic")
A.text(.03,5.55,"thick: the ceiling  $\\psi e^{r s^*}\\!/M(s^*)$",fontsize=9.5,color=ACC)
A.text(.03,5.12,"dashed: measured peak, $t=%.2f$"%tpk,fontsize=9.5,color=INK)
A.text(.03,4.69,"solid: $\\psi$, where the switching wants it",fontsize=9.5,color=BLU)
A.annotate("",xy=(.255,3.75),xytext=(.255,4.45),
           arrowprops=dict(arrowstyle="-|>",color=BLU,lw=1.1))

B=ax[1]
B.axhline(1,color=GRY,lw=.9,ls=":")
B.plot(climb[:,0],climb[:,1],"-",color=INK,lw=1.5,zorder=2)
B.plot(climb[:,0],climb[:,1],"o",ms=5,color=BG,mec=INK,mew=1.4,zorder=3)
i=int(np.argmin(np.abs(climb[:,0]-10.76)))
B.plot(climb[i,0],climb[i,1],"o",ms=8.5,color=ACC,zorder=4)
B.annotate("the paper's own example\nis already at %.0f%% of it"%(100*climb[i,1]),
           (climb[i,0],climb[i,1]),textcoords="offset points",xytext=(12,-46),
           fontsize=9,color=ACC,ha="left",
           arrowprops=dict(arrowstyle="-",color=ACC,lw=.9,
                           connectionstyle="angle3,angleA=0,angleB=80"))
B.set_xscale("log"); B.set_ylim(0,1.15); B.set_xlim(.18,600)
B.set_xlabel("$\\varepsilon=\\bar r_\\psi/\\lambda_1$   (selection speed / mixing speed)",
             color=INK,fontsize=10)
B.set_ylabel("peak imprint  /  ceiling",color=INK,fontsize=10)
B.text(.22,1.10,"it climbs, and never crosses",fontsize=12.5,color=INK,style="italic")

C=ax[2]
C.axhline(0,color=GRY,lw=.9); C.axvline(0,color=GRY,lw=.9)
C.plot(sign[:,0],sign[:,2],"-",color=ACC,lw=4.0,alpha=.42,solid_capstyle="round")
C.plot(sign[:,0],sign[:,1],"o-",ms=4.5,color=INK,lw=1.35,mfc=BG,mew=1.2)
C.set_xlim(-4.2,7.0); C.set_ylim(-.95,1.72)
C.set_xlabel("$\\ln(U^*/U_0)$   —   e-foldings of growth",color=INK,fontsize=10)
C.set_ylabel("shift in mean fitness, at the peak",color=INK,fontsize=10)
C.text(-4.05,1.55,"a shrinking population\nselects backwards",fontsize=12.5,
       color=INK,style="italic",va="top")
C.annotate("enriched for the\nSLOWEST phenotypes",xy=(-2.08,-.598),
           xytext=(-4.05,.30),fontsize=9,color=INK,
           arrowprops=dict(arrowstyle="-|>",color=INK,lw=.9,
                           connectionstyle="angle3,angleA=-70,angleB=10"))
C.annotate("enriched for\nthe fastest",xy=(3.51,1.144),xytext=(3.0,-.72),
           fontsize=9,color=INK,
           arrowprops=dict(arrowstyle="-|>",color=INK,lw=.9,
                           connectionstyle="angle3,angleA=80,angleB=-10"))
C.text(1.45,1.06,"ceiling",color=ACC,fontsize=9.5)
C.text(1.45,.42,"measured",color=INK,fontsize=9.5)

fig.suptitle("Selection's whole budget is the growth.      "
   "$K(s^*)=\\ln(U^*\\!/U_0)$,  with $K$ the cumulant generating function of fitness under $\\psi$",
   color=INK,fontsize=12.5,y=1.02)
fig.text(.995,-.03,"Iris · fire 310 · rebuilt from Fassoni, arXiv:2609.32585",
         ha="right",color=GRY,fontsize=8)
fig.tight_layout()
fig.savefig("selection_budget.png",dpi=155,facecolor=BG,bbox_inches="tight")
print("ok")
