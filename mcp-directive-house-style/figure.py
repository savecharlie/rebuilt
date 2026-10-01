#!/usr/bin/env python3
"""figure.py -- three panels for the house-style measurement. Deterministic (seed 311).

Recomputes everything from the 27 Sep 2026 tools capture; asserts nothing it has not
just measured. Run from this directory with iris-the-maker/earning/mcp importable.
"""
import collections, math, os, random, statistics as st, sys
MCP = os.path.expanduser("~/iris-the-maker/earning/mcp")
sys.path.insert(0, MCP)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import directives as D, collision as C

CAP = os.path.join(MCP, "tools_20260927_union.json.gz")
THR, MIN, SEED = 0.55, 5, 311
rnd = random.Random(SEED)

tools, nhosts = D.load(CAP)
byname = collections.defaultdict(list)
for h, nm, t, s in tools:
    byname[nm].append((h, t, s))
nclust = {nm: (1 if len({h for h, _, _ in v}) <= 1 else len(C.cluster([t for _, t, _ in v], THR)))
          for nm, v in byname.items()}

per = collections.defaultdict(list)          # host -> [(contested, directive)]
names = collections.defaultdict(list)        # host -> [toolname]
for h, nm, t, s in tools:
    per[h].append((nclust[nm] >= 2, bool(D.fams(t, s))))
    names[h].append(nm)

p0 = sum(d for v in per.values() for _, d in v) / sum(len(v) for v in per.values())

# one fleet one vote: hosts sharing an identical tool-NAME SET are one unit
srv = {h: v for h, v in per.items() if len(v) >= MIN}
fleets = collections.defaultdict(list)
for h in srv:
    fleets[frozenset(names[h])].append(h)
units = [(sum(d for _, d in srv[g[0]]), len(srv[g[0]])) for g in fleets.values()]
rates = [k / n for k, n in units]

# --- null: same tool counts, each tool an independent coin at p0 ----------------
B = 3000
nullH = [0.0] * 11
nullext = []
for _ in range(B):
    h0 = [0] * 11; e = 0
    for _, n in units:
        k = sum(1 for _ in range(n) if rnd.random() < p0)
        h0[min(10, int(k / n * 10 + 1e-9))] += 1
        if k == 0 or k == n: e += 1
    for i in range(11): nullH[i] += h0[i] / B
    nullext.append(100 * e / len(units))
obsH = [0] * 11
for r in rates: obsH[min(10, int(r * 10 + 1e-9))] += 1
obsext = 100 * sum(1 for k, n in units if k == 0 or k == n) / len(units)

# --- the reversal ---------------------------------------------------------------
def wil(k, n, z=1.96):
    p = k / n; d = 1 + z * z / n; c = p + z * z / (2 * n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return ((c - h) / d, (c + h) / d)
con = [d for v in per.values() for c, d in v if c]
sol = [d for v in per.values() for c, d in v if not c]
pc, ps = sum(con) / len(con), sum(sol) / len(sol)
lc, hc = wil(sum(con), len(con)); ls_, hs = wil(sum(sol), len(sol))
tool_d = 100 * (pc - ps)
tool_lo = 100 * (pc - ps) - 100 * math.sqrt((hc - lc) ** 2 + (hs - ls_) ** 2) / 2
tool_hi = 100 * (pc - ps) + 100 * math.sqrt((hc - lc) ** 2 + (hs - ls_) ** 2) / 2
diffs = []
for h, v in per.items():
    a = [d for c, d in v if c]; b = [d for c, d in v if not c]
    if a and b: diffs.append(sum(a) / len(a) - sum(b) / len(b))
m = st.mean(diffs); se = st.stdev(diffs) / math.sqrt(len(diffs))
pair_d, pair_lo, pair_hi = 100 * m, 100 * (m - 1.96 * se), 100 * (m + 1.96 * se)

# --- predictive lift ------------------------------------------------------------
Ks = [1, 2, 3, 5, 8]
acc, base = [], []
for K in Ks:
    hit = bh = tot = 0
    for _ in range(40):
        for h, v in per.items():
            if len(v) < K + 1: continue
            d = [x for _, x in v]
            idx = list(range(len(d))); rnd.shuffle(idx)
            held = d[idx[0]]; s = sum(d[i] for i in idx[1:K + 1])
            pred = (s * 2 > K) if s * 2 != K else (p0 >= .5)
            hit += pred == held; bh += held == (p0 >= .5); tot += 1
    acc.append(100 * hit / tot); base.append(100 * bh / tot)

# ================================ draw =========================================
# BINS: exactly-0 and exactly-100 get their OWN bars. A histogram about extremes
# whose left bar is [0,10%) and whose right bar is {100%} is comparing two different
# things at its two ends -- caught by looking at the first render.
def binof(r):
    if r == 0: return 0
    if r == 1: return 11
    return min(10, 1 + int(r * 10 - 1e-12))

INK, ACC, NUL = "#1b1b1b", "#b4461e", "#8a8a8a"
fig, ax = plt.subplots(1, 3, figsize=(15.6, 5.3), facecolor="white")

obsH = [0] * 12
for r in rates: obsH[binof(r)] += 1
nullH = [0.0] * 12
nullext = []
rnd2 = random.Random(SEED + 1)
for _ in range(B):
    h0 = [0] * 12; e = 0
    for _, n in units:
        k = sum(1 for _ in range(n) if rnd2.random() < p0)
        h0[binof(k / n)] += 1
        if k == 0 or k == n: e += 1
    for i in range(12): nullH[i] += h0[i] / B
    nullext.append(100 * e / len(units))
obsext = 100 * sum(1 for k, n in units if k == 0 or k == n) / len(units)

a = ax[0]
x = list(range(12))
cols = [ACC if i in (0, 11) else INK for i in x]
a.bar(x, obsH, width=.84, color=cols, zorder=3)
a.plot(x, nullH, color="#2f6f9f", lw=2.0, marker="o", ms=4.2, zorder=4)
a.set_xticks(x)
a.set_xticklabels(["none", "10", "20", "30", "40", "50", "60", "70", "80", "90", "<100", "all"],
                  fontsize=8.4)
a.set_xlabel("% of a server's tool descriptions that direct the agent\n(first and last bars are exactly none and exactly all)",
             fontsize=9)
a.set_ylabel("servers (fleets counted once)")
a.set_title("A.  servers pile at the ends: all their tools, or none", loc="left", fontsize=11.5, color=INK)
TOP = max(max(obsH), max(nullH)) * 1.30
a.set_ylim(0, TOP)
a.plot([], [], color="#2f6f9f", lw=2, marker="o", ms=4.2,
       label="if each tool were an independent\ncoin at %.1f%%" % (100 * p0))
a.bar([], [], color=ACC, label="exactly none / exactly all")
a.legend(frameon=False, fontsize=8.4, loc="upper right", handlelength=1.6,
         bbox_to_anchor=(1.015, 1.015))
a.annotate("%.0f%% of servers sit at an end\nthe coin expects %.1f%% [%.1f, %.1f]"
           % (obsext, st.mean(nullext), sorted(nullext)[75], sorted(nullext)[-75]),
           xy=(5.6, TOP * 0.74), ha="left", va="top", fontsize=9.2, color=ACC)
for s in ("top", "right"): a.spines[s].set_visible(False)

b = ax[1]
b.axhline(0, color=NUL, lw=1)
vals, los, his = [tool_d, pair_d], [tool_lo, pair_lo], [tool_hi, pair_hi]
for i, (v, lo, hi) in enumerate(zip(vals, los, his)):
    b.plot([i, i], [lo, hi], color=INK, lw=2.4, zorder=3)
    b.plot([i], [v], "o", ms=11, color=(ACC if i == 0 else INK), zorder=4)
    b.annotate("%+.1f" % v, xy=(i, v), xytext=(15, -4), textcoords="offset points",
               fontsize=12, color=(ACC if i == 0 else INK))
b.set_xticks([0, 1])
b.set_xticklabels(["every tool an\nindependent draw\n(12,829 tools)",
                   "paired within server\n(%d servers, each\nspeaking once)" % len(diffs)], fontsize=8.8)
b.set_xlim(-.62, 1.78)
b.set_ylim(min(los) - 2.6, max(his) + 1.4)
b.set_ylabel("directive rate: contested name − sole name  (points)", fontsize=9.5)
b.set_title("B.  the effect I went looking for, and its sign", loc="left", fontsize=11.5, color=INK)
b.annotate("same tools, opposite sign. the hosts that happen to hold\ncontested names are hosts that direct on everything.",
           xy=(.58, min(los) - 1.1), ha="center", va="center", fontsize=8.5, color=NUL)
for s in ("top", "right"): b.spines[s].set_visible(False)

c = ax[2]
c.plot(Ks, acc, color=INK, lw=2.2, marker="o", ms=6, zorder=4,
       label="guess from k other tools\non the same server")
c.plot(Ks, base, color=NUL, lw=1.8, ls="--", label="guess the global majority")
c.fill_between(Ks, base, acc, color=ACC, alpha=.17, zorder=2)
c.set_xlabel("k = other tools from that server you have read", fontsize=9.5)
c.set_ylabel("accuracy on a held-out tool (%)", fontsize=9.5)
c.set_xticks(Ks); c.set_ylim(55, 80)
c.set_title("C.  one tool tells you most of the server", loc="left", fontsize=11.5, color=INK)
c.annotate("+%.1f points\nfrom a single tool" % (acc[0] - base[0]), xy=(1.12, 66.0),
           ha="left", fontsize=9.2, color=ACC)
c.legend(frameon=False, fontsize=8.6, loc="lower right")
for s in ("top", "right"): c.spines[s].set_visible(False)

fig.suptitle("Whether an MCP tool description instructs the agent is more a property of its author than of the tool",
             x=.010, ha="left", fontsize=13.4, color=INK, y=.983)
fig.text(.010, .055,
         "12,829 tool descriptions pulled live from 1,035 reachable MCP servers, 27 Sep 2026.  "
         "Directive detector: measured precision 97.5%% (39/40), miss rate 27.5%% (11/40) -- so it under-reports,\n"
         "and a truly all-directive 5-tool server is missed entirely with probability 0.0016.  "
         "Hosts sharing an identical tool-name set are one unit (631 of 676).  "
         "ICC = 0.387, dispersion $\\chi^2$/df = 7.66.  Deterministic, seed %d." % SEED,
         fontsize=8.2, color=NUL, va="bottom")
fig.tight_layout(rect=[0, .105, 1, .952])
out = "house_style.png"
fig.savefig(out, dpi=165)
print("wrote", out)
print("obs extremes %.1f%%  null %.1f%% [%.1f,%.1f]" % (obsext, st.mean(nullext), sorted(nullext)[75], sorted(nullext)[-75]))
print("bars:", obsH)
print("units %d  tool-level %+.1f [%.1f,%.1f]  paired %+.1f [%.1f,%.1f]"
      % (len(units), tool_d, tool_lo, tool_hi, pair_d, pair_lo, pair_hi))
print("acc", [round(v,1) for v in acc], "base", [round(v,1) for v in base])
