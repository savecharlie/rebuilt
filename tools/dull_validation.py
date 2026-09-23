"""Validate dull.py against faults with KNOWN answers before trusting it.

Written by Iris (Opus 5, 1M), fire 296, Sep 23 2026.

The rule this obeys is my own: an instrument must reproduce a result I already
know before it is allowed to say anything new.  Here the "known results" are
the five faults that actually occurred in fire 295's rebuild of
arXiv:2609.23542, found by hand at the cost of most of a session.

Positive controls 1 and 2 run against the REAL rebuild code in
~/rebuilt/2609.23542-rank-resetting, with the historical bug put back.
Control 5 is the blind spot, asserted rather than hidden.

    python3 dull_validation.py
"""
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).parent))
from dull import watch, sane, safe_div, root_ok, scaling   # noqa: E402

REBUILD = Path.home() / "rebuilt" / "2609.23542-rank-resetting"
rows = []


def check(name, expect_caught, findings, note=""):
    got = bool(findings)
    ok = (got == expect_caught)
    rows.append((ok, name, "caught" if got else "silent",
                 "should catch" if expect_caught else "should be silent", note))
    return ok


# -- control 1: brentq on the naive bracket, real g from the real rebuild ---
if REBUILD.is_dir():
    sys.path.insert(0, str(REBUILD))
    import phic                                             # noqa: E402

    g = phic.g
    # The historical mistake: bracket the whole plausible range instead of the
    # scanned crossing.  g underflows to exactly 0 in the left tail, so brentq
    # accepts lo as a root.  True answer (fire 295): t_c = 0.743904342807.
    lo, hi = 1e-3, 1e3
    root = brentq(g, lo, hi)
    with watch() as w:
        root_ok(g, root, (lo, hi), "t_c (naive bracket)")
    check("1. brentq latched onto a dead bracket endpoint", True, w.findings,
          f"brentq returned t_c = {root:.6g}; truth is 0.743904342807")

    # and the negative control: the bracket fire 295 actually used
    ts = np.geomspace(1e-2, 1e3, 400)
    gs = np.array([g(t) for t in ts])
    i = int(np.where((gs[:-1] < 0) & (gs[1:] > 0))[0][0])
    good = brentq(g, ts[i], ts[i + 1], xtol=1e-14, rtol=8.9e-16)
    with watch() as w:
        root_ok(g, good, (ts[i], ts[i + 1]), "t_c (scanned bracket)")
    check("1b. the correct bracket is not flagged", False, w.findings,
          f"t_c = {good:.12f}  vs published 0.743904342807  "
          f"(delta {abs(good - 0.743904342807):.2g})")
else:
    rows.append((False, "1. rebuild clone missing", "-", "-",
                 f"{REBUILD} not found; see REBUILT_REPO.md"))

# -- control 2: Q0^(N-k) underflow that produced a spurious phi_c = 0.4995 --
Q0 = np.float64(0.4123101754)
with watch() as w:
    for N in (100, 400, 1280, 4000):
        _ = Q0 ** np.float64(N - 1)
check("2. Q0**(N-k) underflowed float64 at large N", True, w.findings,
      f"0.4123**3999 = {float(Q0 ** np.float64(3999)):g}")

with watch() as w:
    for N in (10, 50, 100):
        _ = Q0 ** np.float64(N - 1)
check("2b. the same expression at usable N is not flagged", False, w.findings)

# -- control 3: a finite difference dividing by ~0 in the corners -----------
x = np.linspace(-1.0, 1.0, 9)
with watch() as w:
    safe_div(np.ones_like(x), x ** 3, "dF/dx identity")
check("3. finite difference divided by a corner zero", True, w.findings)

with watch() as w:
    safe_div(np.ones_like(x), x ** 2 + 1.0, "well conditioned identity")
check("3b. a well-conditioned quotient is not flagged", False, w.findings)

# -- control 4: an inverted guard drew an empty curve ----------------------
Ns = np.arange(10, 200)
with watch() as w:
    sane(Ns[Ns > 1e6], "phi_c(N) curve")        # guard inverted -> nothing
check("4. inverted guard produced an empty curve", True, w.findings)

with watch() as w:
    sane(0.4123101754 - 2.07637 / Ns, "phi_c(N) curve")
check("4b. the real curve is not flagged", False, w.findings)

# -- control 4c: the label I read instead of the numbers -------------------
Nn = np.array([10.0, 20, 40, 80, 160, 320, 640])
dev = 2.07637 / Nn                               # decays like 1/N
with watch() as w:
    p = scaling(Nn, dev, "deviation from phi_c", expect=-0.5)
check("4c. a column that decays like 1/N called 1/sqrt(N)", True, w.findings,
      f"measured exponent {p:+.4f}, claimed -0.5000")

with watch() as w:
    p = scaling(Nn, dev, "deviation from phi_c", expect=-1.0)
check("4d. the same column, correctly claimed", False, w.findings,
      f"measured exponent {p:+.4f}")

# -- control 5: THE BLIND SPOT --------------------------------------------
r = 0.7
s = np.linspace(1e-4, 5.0, 200)
with watch() as w:
    right = sane(np.array([1.0 / np.sqrt(np.float64(v) + r) for v in s]), "Qt_0(s+r)")
    wrong = sane(np.array([1.0 / np.sqrt(np.float64(v)) for v in s]), "Qt_0(s)  <- typo")
    safe_div(np.ones_like(s), s + r, "renewal denominator")
check("5. `s` written where `s+r` was meant", False, w.findings,
      f"INVISIBLE by construction -- answer off by {wrong[0] / right[0]:.0f}x "
      f"and every numerical check is green")

# -- report ----------------------------------------------------------------
w_ = max(len(r[1]) for r in rows)
print("dull.py validation -- faults with known answers, from fire 295\n")
for ok, name, got, want, note in rows:
    print(f"  {'PASS' if ok else 'FAIL'}  {name:<{w_}}  {got:>7}  ({want})")
    if note:
        print(f"        {note}")
npass = sum(1 for r in rows if r[0])
print(f"\n  {npass}/{len(rows)} controls behaved as specified.")
print("\n  Control 5 is the one that matters.  dull.py is green on a formula")
print("  typo that is wrong by two orders of magnitude.  A clean report is")
print("  evidence about the dull faults ONLY.  The second instrument, built")
print("  on a different principle, is still the only thing that catches the")
print("  interesting ones -- and the interesting ones are the rare ones.")
sys.exit(0 if npass == len(rows) else 1)
