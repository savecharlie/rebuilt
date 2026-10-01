#!/usr/bin/env python3
"""figure.py -- three panels for the per-token unit. Deterministic at seed 312.

    PYTHONHASHSEED=0 python3 figure.py        -> context_cost.png

A: the unit swap. Same corpus, same hand labels, two denominators.
B: the decoy null. Why the detector's own token gap is not the finding.
C: what survives -- hand-judged directive descriptions are longer where no
   lexicon did the selecting.
"""
from __future__ import annotations
import glob, json, math, random, re, statistics
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import tiktoken
from directives import load, fams

E = tiktoken.get_encoding("o200k_base")
CAP = "tools_20260927_union.json.gz"
INK, DIR, NON, HI = "#1b1b1b", "#b4453c", "#c9c4bb", "#2f6f8f"


def main():
    tools, nsrv = load(CAP)
    tok = {(h, n): len(E.encode((n or "") + "\n" + d)) for h, n, d, _ in tools}
    flag = [t for t in tools if fams(t[2], t[3])]
    T_all = sum(tok.values())
    T_flag = sum(tok[(h, n)] for h, n, _, _ in flag)

    # ---- hand labels
    idx = {(h, n): (d, s) for h, n, d, s in tools}
    rows = []
    for f in sorted(glob.glob("labels_*.jsonl")):
        for line in open(f):
            line = line.strip()
            if not line:
                continue
            r = json.loads(line)
            if "verdict" not in r:
                continue
            k = (r.get("host"), r.get("tool"))
            if k not in idx:
                continue
            d, s = idx[k]
            rows.append({"true": r["verdict"] == "true_positive",
                         "tok": tok[k], "det": bool(fams(d, s))})
    neg = [r for r in rows if not r["det"]]

    # ---- the two-stratum estimates (same arithmetic as tokenweight.py)
    missed = [t for t in tools if not fams(t[2], t[3])]
    T_miss = T_all - T_flag
    s = flag[:]; random.Random(310).shuffle(s); psamp = s[:40]
    plab = [json.loads(l) for l in open("labels_precision_20260930.jsonl")
            if l.strip() and "verdict" in json.loads(l)]
    prec = [(tok[(psamp[r["n"]-1][0], psamp[r["n"]-1][1])],
             r["verdict"] == "true_positive") for r in plab]
    midx = {(t[0], t[1]) for t in missed}
    flo = []
    for l in open("labels_floor_postfix_20260930.jsonl"):
        if not l.strip():
            continue
        r = json.loads(l)
        if "verdict" not in r:
            continue
        k = (r.get("host"), r.get("tool"))
        if k in midx:
            flo.append((tok[k], r["verdict"] == "true_positive"))
    rat = lambda smp: sum(t for t, v in smp if v) / sum(t for t, _ in smp)
    cnt = lambda smp: sum(1 for _, v in smp if v) / len(smp)
    pt = (cnt(prec) * len(flag) + cnt(flo) * len(missed)) / len(tools)
    tw = (rat(prec) * T_flag + rat(flo) * T_miss) / T_all
    rng = random.Random(312)
    bs = []
    for _ in range(20000):
        p = [prec[rng.randrange(len(prec))] for _ in prec]
        f = [flo[rng.randrange(len(flo))] for _ in flo]
        bs.append((rat(p) * T_flag + rat(f) * T_miss) / T_all)
    bs.sort()
    tw_lo, tw_hi = bs[500], bs[19500]

    # ---- decoy null
    WORDS = ["data","return","value","number","list","name","file","text","result",
             "optional","default","format","object","field","string","id","date","type",
             "page","key","api","json","response","request","query","input","output",
             "parameter","code","url","time","user","search","get","set","item","record",
             "status","table","column"]
    texts = [(d, tok[(h, n)]) for h, n, d, _ in tools]
    N = len(texts)
    hits = {w: np.array([bool(re.search(r"\b"+w+r"\b", d, re.I)) for d, _ in texts])
            for w in WORDS}
    toks = np.array([t for _, t in texts])
    target = 100 * len(flag) / N
    real_gap = 100 * T_flag / T_all - target
    r2 = random.Random(312)
    gaps = []
    for _ in range(400):
        sel = np.zeros(N, bool)
        for w in r2.sample(WORDS, len(WORDS)):
            cand = sel | hits[w]
            if 100 * cand.sum() / N > target + 1.5:
                continue
            sel = cand
            if 100 * sel.sum() / N >= target - 1.5:
                break
        p = 100 * sel.sum() / N
        if abs(p - target) > 1.5:
            continue
        gaps.append(100 * toks[sel].sum() / toks.sum() - p)

    # ================================================================ draw
    fig, ax = plt.subplots(1, 3, figsize=(13.4, 4.5))
    plt.rcParams["font.family"] = "DejaVu Sans"

    # -- A: the unit swap
    a = ax[0]
    for i, (lab, share, n_lab) in enumerate(
            [(f"per TOOL\n{len(tools):,} descriptions", pt, f"{100*pt:.1f}%"),
             (f"per TOKEN\n{T_all:,} tokens", tw, f"{100*tw:.1f}%")]):
        a.bar(i, 100, color=NON, edgecolor=INK, linewidth=.8, width=.62)
        a.bar(i, 100 * share, color=DIR, edgecolor=INK, linewidth=.8, width=.62)
        a.text(i - .31, 100 * share + 2.5, n_lab, ha="left", va="bottom",
               fontsize=13, fontweight="bold", color=DIR)
    a.errorbar([0, 1], [100*pt, 100*tw],
               yerr=[[100*pt-44.7, 100*tw-100*tw_lo], [62.5-100*pt, 100*tw_hi-100*tw]],
               fmt="none", ecolor=INK, capsize=5, lw=1.4)
    a.set_xticks([0, 1])
    a.set_xticklabels([f"per TOOL\n{len(tools):,} descriptions",
                       f"per TOKEN\n{T_all:,} tokens"], fontsize=9)
    a.set_ylim(0, 112)
    a.set_ylabel("share that instructs the agent (hand labels)", fontsize=9)
    a.set_title("A  one corpus, two denominators", fontsize=10, loc="left", pad=8)
    a.annotate("", xy=(1, 100*tw), xytext=(0, 100*pt),
               arrowprops=dict(arrowstyle="->", color=HI, lw=1.6,
                               connectionstyle="arc3,rad=-0.25"))
    a.text(.5, 86, f"+{100*(tw-pt):.0f} points", ha="center", fontsize=10,
           color=HI, fontweight="bold")
    for sp in ("top", "right"):
        a.spines[sp].set_visible(False)

    # -- B: the decoy null
    b = ax[1]
    b.hist(gaps, bins=np.arange(6, 18.01, .5), color=NON, edgecolor=INK, linewidth=.6)
    b.axvline(real_gap, color=DIR, lw=2.4)
    b.text(real_gap + .35, b.get_ylim()[1]*.94, f"real detector\n+{real_gap:.1f} pts",
           ha="left", va="top", fontsize=9, color=DIR, fontweight="bold")
    b.text(statistics.mean(gaps), b.get_ylim()[1]*.58,
           f"400 decoy lexicons\nmean +{statistics.mean(gaps):.1f}  sd {statistics.stdev(gaps):.1f}",
           ha="center", va="center", fontsize=9, color=INK,
           bbox=dict(fc="white", ec=INK, lw=.6, alpha=.92, pad=3))
    b.set_xlabel("token-share minus tool-share, percentage points", fontsize=9)
    b.set_ylabel("decoy lexicons", fontsize=9)
    b.set_title("B  the gap a word list gets for free", fontsize=10, loc="left", pad=8)
    b.set_xlim(6, 19.5)
    for sp in ("top", "right"):
        b.spines[sp].set_visible(False)

    # -- C: what survives
    c = ax[2]
    groups = [("detector says NO\n(no lexicon selected these)",
               [r["tok"] for r in neg if r["true"]],
               [r["tok"] for r in neg if not r["true"]]),
              ("all 266 labelled rows",
               [r["tok"] for r in rows if r["true"]],
               [r["tok"] for r in rows if not r["true"]])]
    rs = np.random.default_rng(312)
    for gi, (lab, t_, f_) in enumerate(groups):
        for j, (vals, col) in enumerate(((f_, NON), (t_, DIR))):
            x = gi * 2 + j + rs.normal(0, .07, len(vals))
            c.scatter(x, vals, s=11, color=col, edgecolor=INK, linewidth=.25,
                      alpha=.75, zorder=3)
            m = statistics.median(vals)
            c.plot([gi*2+j-.3, gi*2+j+.3], [m, m], color=INK, lw=2.2, zorder=4)
            c.text(gi*2+j, 700, f"med {m:.0f}", ha="center", fontsize=8, color=INK)
    c.set_yscale("log")
    c.set_ylim(3, 2600)
    c.set_xticks([0, 1, 2, 3])
    c.set_xticklabels(["not\ndirective", "directive", "not\ndirective", "directive"],
                      fontsize=8)
    c.set_ylabel("description length, o200k tokens (log)", fontsize=9)
    c.set_title("C  what survives: judged by hand", fontsize=10, loc="left", pad=8)
    c.text(0.5, 1500, groups[0][0], ha="center", va="center", fontsize=8, color=HI)
    c.text(2.5, 1500, groups[1][0], ha="center", va="center", fontsize=8, color=HI)
    c.axvline(1.5, color=INK, lw=.6, ls=":", alpha=.5)
    neg_t = [r["tok"] for r in neg if r["true"]]
    neg_f = [r["tok"] for r in neg if not r["true"]]
    c.text(0.5, 4.6, f"{statistics.mean(neg_t)/statistics.mean(neg_f):.2f}x   z=+2.93",
           ha="center", fontsize=9, color=DIR, fontweight="bold")
    allt = [r["tok"] for r in rows if r["true"]]; allf = [r["tok"] for r in rows if not r["true"]]
    c.text(2.5, 4.6, f"{statistics.mean(allt)/statistics.mean(allf):.2f}x   z=+7.75",
           ha="center", fontsize=9, color=DIR, fontweight="bold")
    for sp in ("top", "right"):
        c.spines[sp].set_visible(False)

    fig.suptitle("What an MCP tool catalogue spends of an agent's context window  "
                 "— 1,035 live servers, 12,829 tool descriptions, 27 Sep 2026",
                 fontsize=10.5, y=1.00, x=.5)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig("context_cost.png", dpi=170)
    print("wrote context_cost.png")
    print(f"  A: per-tool {100*pt:.1f}%  per-token {100*tw:.1f}% [{100*tw_lo:.1f}, {100*tw_hi:.1f}]")
    print(f"  B: real {real_gap:+.1f}  decoys n={len(gaps)} mean {statistics.mean(gaps):+.1f} sd {statistics.stdev(gaps):.1f}")
    print(f"  C: neg stratum {len(neg)} rows, ratio {statistics.mean(neg_t)/statistics.mean(neg_f):.2f}x")


if __name__ == "__main__":
    main()
