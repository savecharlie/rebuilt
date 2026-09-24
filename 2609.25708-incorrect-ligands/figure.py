import numpy as np, math, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ising import Cluster
cl=Cluster(5); pin=cl.pin_centre_up()
J=0.55; D=4.0
h_of_c=lambda c:(math.log(c)+1)/2
h0=h_of_c(0.32)
fig,(ax,bx)=plt.subplots(1,2,figsize=(12.6,4.9))

# ---- A: the claim is "activation is a bigger field than the whole concentration axis"
hs=np.linspace(-1.2,2.6,420)
ax.plot(hs,[cl.mean_bound(J,h) for h in hs],lw=2.4,color="#1b3a6b",label="no correct ligand")
ax.plot(hs,[pin.mean_bound(J,h) for h in hs],lw=2.0,color="#c0392b",ls="--",
        label="one correct ligand (centre)")
lo,hi=h_of_c(0.15),h_of_c(0.50)
ax.axvspan(lo,hi,color="#f0c419",alpha=.30,zorder=0)
ax.text((lo+hi)/2,0.72,"their ENTIRE\nc axis\n0.15–0.50",ha="center",va="center",fontsize=8.5,color="#7a5c00")
for a,col in ((0.0,"#1b3a6b"),(0.5,"#2e8b57"),(1.0,"#6a0dad")):
    h=h0+D*a/2; y=cl.mean_bound(J,h)
    ax.plot([h],[y],"o",ms=8,color=col,zorder=5)
    ax.annotate(f"a={'0' if a==0 else ('1/2' if a==.5 else '1')}\n{y:.3f}",(h,y),
                textcoords="offset points",xytext=(8,-4 if a==0 else -22),fontsize=9,color=col)
ax.annotate("", xy=(h0+2.0,0.5), xytext=(h0,0.5),
            arrowprops=dict(arrowstyle="-|>",lw=2.2,color="#6a0dad"))
ax.text(h0+1.0,0.44,"one activation step = Δε/2 = 1.0\nfull activation = 2.0",
        ha="center",va="top",fontsize=9,color="#6a0dad")
ax.set_xlabel("effective field  h = (ln c + 1)/2"); ax.set_ylabel("equilibrium bound fraction")
ax.set_title(f"A.  exact 5×5 equilibrium, $J_b$={J}\nthe ratchet dwarfs the concentration axis",fontsize=10.5)
ax.set_ylim(-0.03,1.16); ax.legend(fontsize=8.5,loc="upper left",framealpha=.92); ax.grid(alpha=.25)

# ---- B: the free-energy climb, and what the seed is worth
n=np.arange(26)
Fu=cl.Fn(J,h0); Fs=pin.Fn(J,h0)
Fu=Fu-np.nanmin(Fu); Fs=Fs-np.nanmin(Fs)
bx.plot(n,Fu,"o-",ms=4,lw=1.9,color="#1b3a6b",label="no correct ligand")
bx.plot(n,Fs,"s--",ms=4,lw=1.9,color="#c0392b",label="one correct ligand")
bx.axvline(10,color="k",lw=.9,ls=":"); bx.text(12.5,1.2,"their $N_a$=10\nthreshold",fontsize=8.5)
bx.annotate("",xy=(10,Fs[10]),xytext=(10,Fu[10]),arrowprops=dict(arrowstyle="<|-|>",lw=1.8,color="#2e8b57"))
bx.text(11.2,(Fu[10]+Fs[10])/2-2.6,f"  {Fu[10]-Fs[10]:.2f} $k_BT$\n  = {(Fu[10]-Fs[10])/J:.1f}·$J_b$\n  (rate ×{math.exp(Fu[10]-Fs[10]):.0f})",
        fontsize=9,color="#2e8b57",ha="left",va="center")
bx.set_xlabel("n  (bound receptors)"); bx.set_ylabel("F(n) − F(min)   [$k_BT$]")
bx.set_title("B.  the exact climb at their operating point (c=0.32)\nno second minimum: at a=0 there is nothing to nucleate INTO",fontsize=10.5)
bx.legend(fontsize=8.5); bx.grid(alpha=.25); bx.set_xlim(-0.5,25.5)
plt.tight_layout(); plt.savefig("ratchet.png",dpi=155)
print("wrote ratchet.png")
print(f"  seed advantage at n=10: {Fu[10]-Fs[10]:.4f} kT = {(Fu[10]-Fs[10])/J:.3f} J_b, rate x{math.exp(Fu[10]-Fs[10]):.1f}")
