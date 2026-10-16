#!/usr/bin/env python3
"""fig6.py -- the trapped state is three walls wearing a cloud of short-lived pairs.

WHAT THIS CORRECTS AND WHAT IT DOES NOT. NOTES.md says survivors are 3 or 6,
always, never 4 or 5, because a same-charge set needs n = 0 mod 3. That is a
claim about the PERSISTENT wall set and it stands. What a window-free reading
adds is that the persistent set is dressed: charge-1 walls continually bud
charge-2 pairs that reabsorb within a few steps. Mean episode 4.9 steps, median
3, max 31, 224 episodes in 1,600 steps on one ring. A 15-block smoothing window
reports three walls in 98% of frames because the excitations are both narrow and
brief, so the old number is what that instrument can see and not an error in it.

The same cloud is absent on the other chirality: (2,2,2) rings reach four walls
27% of the time and never five. The split channel 1 -> 2+2 is open and 2 -> 1+1
is shut, which is the inversion of registered P19/P21.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

from split import plateaus
from tri import Kit
from wallcount import locked_ring as wl
from wallcount import step as ws
from wallcount import walls as ww

kit = Kit((36, 37, 5))
N, T0, T, W = 300, 2_000, 1_200, 15
CM = ListedColormap(["#d9d4cc", "#2d4a6b", "#b5462f", "#d9a326"])


def ring(beats, seed):
    rng = np.random.default_rng(seed)
    s = wl(N, rng)
    for i, off in enumerate(beats):
        b = wl(N, rng)
        for _ in range(off):
            b = ws(b)
        lo, hi = i * N // 3, (i + 1) * N // 3
        s[3 * lo:3 * hi] = b[3 * lo:3 * hi]
    return s


def collect(beats, seed):
    s = ring(beats, seed)
    for _ in range(T0):
        s = kit.step(s)
    field, npl, nw, short = [], [], [], []
    for t in range(T0, T0 + T):
        field.append(kit.labels(s, t).astype(np.int16))
        pl = plateaus(kit, s, t)
        npl.append(len(pl) if len(pl) >= 2 else 0)
        nw.append(int(ww(s, t, W)[1]))
        short.append(min(pl, key=lambda r: r[1])[2] if len(pl) >= 4 else -1)
        s = kit.step(s)
    return np.array(field), np.array(npl), np.array(nw), np.array(short)


def episodes(npl, thresh=4):
    out, cur = [], 0
    for n in npl:
        if n >= thresh:
            cur += 1
        elif cur:
            out.append(cur)
            cur = 0
    return np.array(out)


f1, npl1, nw1, sh1 = collect((0, 1, 2), 7003)
f2, npl2, nw2, _ = collect((0, 2, 1), 7003)
ep = episodes(npl1)

fig = plt.figure(figsize=(16, 8.4))
gs = fig.add_gridspec(2, 3, width_ratios=[1.25, 0.85, 1.0], height_ratios=[1, 1],
                      hspace=0.34, wspace=0.26)

# (a) the whole ring
ax = fig.add_subplot(gs[:, 0])
ax.imshow(f1 + 1, aspect="auto", interpolation="nearest", cmap=CM, vmin=0, vmax=3,
          origin="lower", extent=[0, N, T0, T0 + T])
ax.set_xlabel("block around the ring")
ax.set_ylabel("step  t")
ax.set_title("(a) the trapped state, charges (1,1,1)\nraw phase labels, no window",
             fontsize=11)
ax.text(0.02, 0.985, "grey = no beat-exclusive form\nblue/red/gold = phase 0/1/2",
        transform=ax.transAxes, va="top", fontsize=8,
        bbox=dict(fc="white", ec="0.7", alpha=0.93, pad=3))

# (b) zoom on the single clearest bud: the frame with the most walls, and the
# short plateau's own position, so the zoom is aimed by the data not by hand.
j = int(np.argmax(npl1))
cx = int(sh1[j]) if sh1[j] >= 0 else N // 2
x0, x1 = max(0, cx - 14), min(N, cx + 14)
t0, t1 = max(0, j - 14), min(T, j + 14)
bx = fig.add_subplot(gs[0, 1])
bx.imshow(f1[t0:t1, x0:x1] + 1, aspect="auto", interpolation="nearest", cmap=CM,
          vmin=0, vmax=3, origin="lower", extent=[x0, x1, T0 + t0, T0 + t1])
bx.axhline(T0 + j, color="k", lw=0.9, ls=":")
bx.annotate("a strip of the third phase,\nborn and reabsorbed",
            xy=(cx + 0.5, T0 + j), xytext=(0.04, 0.06), textcoords="axes fraction",
            fontsize=8, color="k",
            arrowprops=dict(arrowstyle="->", lw=1.1, color="k"),
            bbox=dict(fc="white", ec="0.7", alpha=0.95, pad=2.5))
bx.set_title(f"(b) one bud, {x1-x0} blocks x {t1-t0} steps\n"
             f"{npl1[j]} walls at t={T0+j}, where w={W} reads {nw1[j]}", fontsize=11)
bx.set_xlabel("block")
bx.set_ylabel("step  t")

# (c) episode lifetimes
cx_ = fig.add_subplot(gs[1, 1])
u, c = np.unique(ep, return_counts=True)
cx_.bar(u, c, color="#2d4a6b", width=0.85)
cx_.set_yscale("log")
cx_.set_xlabel("episode length (steps with >=4 walls)")
cx_.set_ylabel("count")
cx_.set_title(f"(c) the buds are brief\n{len(ep)} episodes, mean {ep.mean():.1f}, "
              f"median {np.median(ep):.0f}, max {ep.max()} steps", fontsize=11)

# (d) the two counts, both chiralities
dx = fig.add_subplot(gs[:, 2])
bins = np.arange(-0.5, 9, 1)
for arr, lab, col, off in ((npl1, "(1,1,1) raw labels", "#2d4a6b", -0.3),
                           (nw1, f"(1,1,1) smoothed w={W}", "#b5462f", -0.1),
                           (npl2, "(2,2,2) raw labels", "#d9a326", 0.1),
                           (nw2, f"(2,2,2) smoothed w={W}", "#7a7a7a", 0.3)):
    h, _ = np.histogram(arr, bins=bins)
    dx.barh(np.arange(9) + off, 100 * h / h.sum(), height=0.2, color=col, label=lab)
dx.set_yticks(range(9))
dx.set_ylim(1.5, 8.5)
dx.set_ylabel("walls detected")
dx.set_xlabel("% of frames")
dx.legend(fontsize=8.5, loc="lower right", framealpha=0.95)
dx.grid(axis="x", alpha=0.3)
dx.set_title("(d) what each instrument sees\nthe cloud is on the (1,1,1) chirality only",
             fontsize=11)

fig.suptitle("arXiv:2610.09879, triple (36,37,5) — the trapped state is three walls "
             "dressed in short-lived pairs (fire 344)", fontsize=13)
fig.savefig("section6_window.png", dpi=120, bbox_inches="tight")
print(f"(1,1,1) raw {dict(zip(*[x.tolist() for x in np.unique(npl1, return_counts=True)]))}")
print(f"(1,1,1) w{W} {dict(zip(*[x.tolist() for x in np.unique(nw1, return_counts=True)]))}")
print(f"(2,2,2) raw {dict(zip(*[x.tolist() for x in np.unique(npl2, return_counts=True)]))}")
print(f"(2,2,2) w{W} {dict(zip(*[x.tolist() for x in np.unique(nw2, return_counts=True)]))}")
print(f"episodes n={len(ep)} mean {ep.mean():.2f} median {np.median(ep):.0f} max {ep.max()}")
print("wrote section6_window.png")
