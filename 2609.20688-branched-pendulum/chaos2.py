import numpy as np, sys
from chaos import lyap
th = float(sys.argv[1])
print(th, " ".join(f"T={T}:{lyap(th, T=T):+.3f}" for T in (60, 240)), flush=True)
