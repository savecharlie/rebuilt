"""Do the three traditional AMOC fingerprints agree with EACH OTHER in the
observations?  They are all claimed to measure one quantity, so this test
needs no model at all: if they disagree, at least two of them are wrong.
"""
import numpy as np, ersst, hadisst, fingerprints as fp

def rolling(x, n):
    if n == 1: return x.astype(float)
    k = np.ones(n)/n
    out = np.convolve(x.astype(float), k, 'valid')
    return out

def report(name, cls):
    e, f = fp.build(cls())
    # put all three on a common annual axis (DP is calendar-annual, SG winter-year)
    ys, ya = f['years'], f['years_a']
    lo, hi = max(ys[0], ya[0], 1871), min(ys[-1], ya[-1], 2024)
    ms, ma = (ys>=lo)&(ys<=hi), (ya>=lo)&(ya<=hi)
    S = {'SG-G': np.asarray(f['SG-G'][ms]),
         'SG-NH': np.asarray(f['SG-NH'][ms]),
         'DP':   np.asarray(f['DP'][ma])}
    n = min(len(v) for v in S.values()); S = {k:v[:n] for k,v in S.items()}
    print(f"\n===== {name}  {lo}-{hi}, n={n} =====")
    for smooth in (1, 5, 30):
        keys = list(S)
        print(f"  {smooth:2d}-yr mean   " + "  ".join(f"{a}~{b}" for i,a in enumerate(keys) for b in keys[i+1:]))
        for label, prep in (('raw     ', lambda v: v),
                            ('detrended', lambda v: v - np.polyval(np.polyfit(np.arange(len(v)), v, 1), np.arange(len(v))))):
            vals = {k: rolling(prep(v), smooth) for k, v in S.items()}
            rs = [np.corrcoef(vals[a], vals[b])[0,1] for i,a in enumerate(keys) for b in keys[i+1:]]
            print(f"     {label}  " + "   ".join(f"{r:+6.3f}" for r in rs))

for nm, c in [('ERSSTv5', ersst.ERSST), ('HadISST', hadisst.HadISST)]:
    report(nm, c)
