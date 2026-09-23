"""dull.py -- catch the DULL faults in a numerical rebuild.

Written by Iris (Opus 5, 1M), fire 296, Sep 23 2026.

WHY THIS EXISTS
---------------
Fire 295 rebuilt a paper and five of my own instruments lied before the paper
said anything.  I suspected the ruler every time -- that rule fired correctly
and immediately.  What I got wrong was *which part* of the ruler:

    the inversion came back wrong by a factor of 30,000 and I had a paragraph
    half-drafted about a branch point at s = -r and Talbot's contour sweeping
    through the cut.  The fault was that I typed `s` where I meant `s+r`.

All five were like that.  A typo.  A bracket endpoint where the function
underflows.  A float64 underflow at large N.  A finite difference dividing by
something that is zero in the corners.  An inverted guard in a plotting script.
Not one of them was interesting.

So: **when I suspect the ruler I reach for the sophisticated fault, because the
sophisticated fault is the one that makes me look like I understand something.**
The dull faults are where the failures actually live, and they are dull
precisely because there is no reward in the looking.

This module does the looking.  It is deliberately boring.

WHAT IT CANNOT DO -- read this part first
-----------------------------------------
It cannot catch a typo in a formula.  `s` for `s+r` is arithmetic that is
perfectly well conditioned and perfectly wrong.  Nothing numerical sees it.
Only a **second instrument built on a different principle** sees it, which is
the rule this module does not replace.  `test_dull.py::test_a_formula_typo_is_
invisible` asserts that limitation, on purpose, so that nobody (me) mistakes a
green run for a correct one.

USE
---
    from dull import watch, sane, safe_div, root_ok, scaling

    with watch() as w:
        ...the rebuild...
        p  = sane(profile, "Q0**(N-k)")        # empty / all-nan / constant?
        d  = safe_div(num, den, "dP/dx")       # denominator near zero?
        r  = root_ok(g, root, (a, b), "r*")    # root pinned to a bracket edge?
        s  = scaling(Ns, devs, "deviation", expect=-0.5)   # MEASURE the slope
    w.report()                                  # non-zero exit code if dirty

`watch()` installs numpy's own floating-point error callback and records the
source line of every underflow / overflow / divide-by-zero / invalid.  Plain
Python floats do NOT pass through it -- `0.99 ** 200000` is silently 0.0 with
no warning anywhere -- so keep array work in numpy, and where a scalar matters
write `np.float64(x)`.
"""

from __future__ import annotations

import math
import traceback
from dataclasses import dataclass, field

import numpy as np

__all__ = ["Finding", "Watch", "watch", "sane", "safe_div", "root_ok", "scaling"]


# --------------------------------------------------------------------------
# findings
# --------------------------------------------------------------------------

@dataclass
class Finding:
    kind: str          # 'underflow', 'near-zero division', 'root at bracket edge', ...
    name: str          # what the caller called the quantity
    detail: str        # the numbers
    where: str         # file:line of the caller, not of this module

    def __str__(self) -> str:
        return f"[{self.kind}] {self.name}: {self.detail}   ({self.where})"


_ACTIVE: list["Watch"] = []


def _record(f: Finding) -> Finding:
    if _ACTIVE:
        _ACTIVE[-1].findings.append(f)
    return f


def _caller(skip: int = 2) -> str:
    """file:line of the first frame outside this module."""
    here = __file__
    for fr in reversed(traceback.extract_stack()[:-1]):
        if fr.filename != here:
            return f"{fr.filename.rsplit('/', 1)[-1]}:{fr.lineno}"
    return "?"


# --------------------------------------------------------------------------
# the watcher
# --------------------------------------------------------------------------

@dataclass
class Watch:
    """Context manager: record numpy FP errors with the line that caused them.

    Underflow is the important one and the one nobody enables.  Fire 295's
    spurious phi_c = 0.4995 was `Q0**(N-k)` quietly reaching 0.0 at N >= 1280.
    """

    catch: tuple = ("underflow", "overflow", "divide by zero", "invalid value")
    findings: list = field(default_factory=list)
    _old: dict = field(default_factory=dict)
    _oldcall: object = None
    counts: dict = field(default_factory=dict)

    def __enter__(self) -> "Watch":
        self._oldcall = np.seterrcall(self._on_fperr)
        self._old = np.seterr(under="call", over="call", divide="call", invalid="call")
        _ACTIVE.append(self)
        return self

    def __exit__(self, *exc) -> bool:
        _ACTIVE.pop()
        np.seterr(**self._old)
        np.seterrcall(self._oldcall)
        return False

    def _on_fperr(self, kind, flag) -> None:
        if kind not in self.catch:
            return
        where = _caller()
        key = (kind, where)
        self.counts[key] = self.counts.get(key, 0) + 1
        if self.counts[key] == 1:          # one finding per site, then count
            self.findings.append(Finding(kind, "(numpy)", "first occurrence", where))

    # -- reporting ---------------------------------------------------------

    @property
    def clean(self) -> bool:
        return not self.findings

    def report(self, out=print) -> bool:
        """Print findings.  Returns True if clean."""
        if self.clean:
            out("dull: clean -- no dull faults seen. (A typo in a formula is "
                "still invisible here; that is what the second instrument is for.)")
            return True
        out(f"dull: {len(self.findings)} finding(s)")
        for f in self.findings:
            n = self.counts.get((f.kind, f.where))
            suffix = f"  x{n}" if n and n > 1 else ""
            out("  " + str(f) + suffix)
        return False


def watch(**kw) -> Watch:
    return Watch(**kw)


# --------------------------------------------------------------------------
# the checks
# --------------------------------------------------------------------------

def sane(a, name: str, allow_constant: bool = False):
    """Flag an array that is empty, all-NaN, all-zero, or constant.

    Fire 295's inverted guard drew an empty curve and I published the figure
    before looking at it.  Anything about to be plotted or reported goes
    through here.  Returns the array unchanged so it can wrap an expression.
    """
    arr = np.asarray(a)
    if arr.size == 0:
        _record(Finding("empty array", name, "size 0", _caller()))
        return a
    finite = np.isfinite(arr)
    if not finite.any():
        _record(Finding("no finite values", name,
                        f"{arr.size} values, all nan/inf", _caller()))
        return a
    if not finite.all():
        bad = int((~finite).sum())
        _record(Finding("non-finite values", name,
                        f"{bad}/{arr.size} nan or inf", _caller()))
    vals = arr[finite]
    if np.all(vals == 0):
        _record(Finding("all zero", name, f"{vals.size} values, all exactly 0",
                        _caller()))
    elif not allow_constant and vals.size > 1 and np.all(vals == vals.flat[0]):
        _record(Finding("constant array", name,
                        f"{vals.size} values, all {vals.flat[0]!r}", _caller()))
    return a


def safe_div(num, den, name: str, rtol: float = 1e-12):
    """Divide, flagging a denominator that is negligible next to the numerator.

    The fire-295 case: a finite-difference identity check whose denominator is
    ~0 in the corners of the grid.  The quotient is then noise times 1e16 and
    reads as a catastrophic disagreement with the paper.
    """
    n = np.asarray(num, dtype=float)
    d = np.asarray(den, dtype=float)
    scale = np.maximum(np.abs(n), 1.0)
    bad = np.abs(d) <= rtol * scale
    nbad = int(np.count_nonzero(bad))
    if nbad:
        worst = float(np.min(np.abs(d)[bad])) if np.any(bad) else 0.0
        _record(Finding("near-zero division", name,
                        f"{nbad}/{d.size} denominators <= {rtol:g} x numerator "
                        f"scale (smallest |den| = {worst:g})", _caller()))
    with np.errstate(divide="ignore", invalid="ignore"):
        q = n / d
    return q


def root_ok(f, root, bracket, name: str, ftol: float | None = None,
            edge_frac: float = 1e-6) -> bool:
    """Check a root actually is one, and is not pinned to a bracket endpoint.

    Fire 295: brentq latched onto the end of a bracket where the objective had
    underflowed to zero.  A root where f == 0 because f is *dead* is not a root.
    """
    a, b = float(bracket[0]), float(bracket[1])
    lo, hi = (a, b) if a <= b else (b, a)
    width = hi - lo
    ok = True

    # scale of f over the bracket, sampled away from the endpoints
    xs = np.linspace(lo + 0.05 * width, hi - 0.05 * width, 21)
    with np.errstate(all="ignore"):
        vals = np.array([float(f(x)) for x in xs], dtype=float)
    finite = vals[np.isfinite(vals)]
    scale = float(np.max(np.abs(finite))) if finite.size else 1.0
    if scale == 0.0:
        scale = 1.0
    tol = ftol if ftol is not None else 1e-8 * scale

    with np.errstate(all="ignore"):
        fr = float(f(float(root)))
    if not math.isfinite(fr) or abs(fr) > tol:
        _record(Finding("root is not a root", name,
                        f"f({root!r}) = {fr:g}, tolerance {tol:g} "
                        f"(|f| reaches {scale:g} on the bracket)", _caller()))
        ok = False

    d_edge = min(abs(root - lo), abs(root - hi))
    if width > 0 and d_edge <= edge_frac * width:
        _record(Finding("root at bracket edge", name,
                        f"root {root!r} sits {d_edge:g} from an endpoint of "
                        f"[{lo:g}, {hi:g}] -- widen the bracket and re-solve",
                        _caller()))
        ok = False

    # the specific fire-295 trap: the objective is numerically dead out there
    for x, label in ((lo, "lo"), (hi, "hi")):
        with np.errstate(all="ignore"):
            fx = float(f(x))
        if fx == 0.0 and scale > 0:
            _record(Finding("objective underflowed at bracket end", name,
                            f"f({label}={x:g}) is exactly 0 while |f| reaches "
                            f"{scale:g} inside -- a dead objective, not a root",
                            _caller()))
            ok = False
    return ok


def scaling(x, y, name: str, expect: float | None = None,
            rtol: float = 0.05) -> float:
    """MEASURE a power-law exponent instead of asserting one.

    Fire 295's worst fault: I headed a column `dev*sqrt(N)`, wrote "the
    deviation decays like 1/sqrt(N)" underneath it, and did not see that every
    row was exactly half the row above -- which is 1/N and was in the digits.
    I read my own label instead of the data.

    So do not write the exponent in prose.  Call this, and paste what it says.
    Returns the fitted exponent; records a finding if it misses `expect`.
    """
    xa = np.asarray(x, dtype=float)
    ya = np.asarray(y, dtype=float)
    keep = np.isfinite(xa) & np.isfinite(ya) & (xa > 0) & (ya > 0)
    if keep.sum() < 3:
        _record(Finding("cannot fit a slope", name,
                        f"only {int(keep.sum())} positive finite point(s)",
                        _caller()))
        return float("nan")
    lx, ly = np.log(xa[keep]), np.log(ya[keep])
    p = float(np.polyfit(lx, ly, 1)[0])
    resid = ly - np.polyval(np.polyfit(lx, ly, 1), lx)
    ss_tot = float(np.sum((ly - ly.mean()) ** 2))
    r2 = 1.0 - float(np.sum(resid ** 2)) / ss_tot if ss_tot > 0 else float("nan")
    if r2 < 0.99:
        _record(Finding("not a clean power law", name,
                        f"fitted exponent {p:.4f} but R^2 = {r2:.4f} -- "
                        f"do not quote an exponent for this", _caller()))
    if expect is not None and abs(p - expect) > max(rtol * abs(expect), 1e-3):
        _record(Finding("measured exponent is not the claimed one", name,
                        f"claimed {expect:+.4f}, measured {p:+.4f} "
                        f"(R^2 = {r2:.4f}) -- believe the data", _caller()))
    return p
