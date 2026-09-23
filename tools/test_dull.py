"""Tests for dull.py -- every case is a real fault from fire 295's rebuild.

Written by Iris (Opus 5, 1M), Sep 23 2026.

Five instruments lied that day and I reached for the sophisticated explanation
every single time.  Each test below reconstructs one of them.  The last test
asserts what this module CANNOT see, which is the most important one.
"""

import math
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).parent))
from dull import watch, sane, safe_div, root_ok, scaling   # noqa: E402


def kinds(w):
    return [f.kind for f in w.findings]


# -- fault 3: Q0**(N-k) underflowed float64 at N >= 1280 --------------------

def test_underflow_in_a_power_is_caught_with_its_line():
    with watch() as w:
        Q0 = np.float64(0.98)
        val = Q0 ** np.float64(40000)      # -> 0.0, silently, in plain python
    assert val == 0.0
    assert "underflow" in kinds(w)
    assert w.findings[0].where.startswith("test_dull.py:")
    assert not w.clean


def test_plain_python_floats_are_invisible_and_the_docstring_says_so():
    """Honest limit: the callback is numpy's, so `0.98 ** 40000` slips past."""
    with watch() as w:
        val = 0.98 ** 40000
    assert val == 0.0
    assert w.clean                          # nothing seen -- by construction
    assert "np.float64" in __import__("dull").__doc__


# -- fault 4: a finite-difference identity dividing by ~0 in the corners ----

def test_near_zero_denominator_is_flagged():
    x = np.linspace(-1, 1, 5)
    num = np.ones_like(x)
    den = x ** 3                            # exactly 0 at the middle sample
    with watch() as w:
        q = safe_div(num, den, "dP/dx identity")
    assert "near-zero division" in kinds(w)
    assert not np.isfinite(q[2])
    assert "1/5" in w.findings[0].detail


def test_ordinary_division_is_silent():
    with watch() as w:
        q = safe_div(np.array([2.0, 4.0]), np.array([1.0, 2.0]), "ratio")
    assert w.clean
    assert np.allclose(q, [2.0, 2.0])


# -- fault 2: brentq latched onto a bracket end where the objective is dead --

def test_root_pinned_to_a_dead_bracket_endpoint():
    # The fire-295 shape: an objective built from powers that both underflow
    # to exactly 0 far out.  g is strictly POSITIVE wherever it is alive -- it
    # has no root at all -- but brentq is handed g(hi) == 0.0 and returns hi.
    def g(x):
        with np.errstate(under="ignore"):
            return float(np.float64(0.5) ** np.float64(x)
                         - np.float64(0.4) ** np.float64(x))

    lo, hi = 1.0, 1200.0
    root = brentq(g, lo, hi)
    assert root == hi                       # brentq handed back the endpoint
    with watch() as w:
        ok = root_ok(g, root, (lo, hi), "r*")
    assert not ok
    assert "objective underflowed at bracket end" in kinds(w)
    assert "root at bracket edge" in kinds(w)


def test_a_genuine_root_passes_quietly():
    f = lambda x: math.cos(x) - x
    root = brentq(f, 0.0, 1.5)
    with watch() as w:
        ok = root_ok(f, root, (0.0, 1.5), "cos x = x")
    assert ok
    assert w.clean
    assert abs(root - 0.7390851332151607) < 1e-12


def test_a_root_at_the_very_edge_is_flagged_even_if_f_is_alive():
    f = lambda x: x - 1.0
    with watch() as w:
        ok = root_ok(f, 1.0, (1.0, 5.0), "edge root")
    assert not ok
    assert "root at bracket edge" in kinds(w)


# -- fault 5: an inverted guard drew an empty curve -------------------------

def test_empty_and_degenerate_arrays_are_flagged_before_plotting():
    with watch() as w:
        sane(np.array([]), "curve")
        sane(np.array([np.nan, np.nan]), "phi_c(N)")
        sane(np.zeros(7), "residual")
        sane(np.full(5, 3.3), "threshold")
    assert kinds(w) == ["empty array", "no finite values", "all zero",
                        "constant array"]


def test_a_real_curve_passes():
    with watch() as w:
        sane(np.linspace(0.3, 0.5, 40), "phi_c(N)")
    assert w.clean


def test_sane_returns_its_argument_so_it_can_wrap_an_expression():
    a = np.arange(5.0)
    with watch():
        assert sane(a, "x") is a


# -- the worst one: I read my own column label instead of the numbers -------

def test_scaling_contradicts_a_wrong_claim():
    """The dev*sqrt(N) column.  Every row was half the one above -- 1/N --
    and I wrote '1/sqrt(N)' under it."""
    N = np.array([10.0, 20, 40, 80, 160, 320])
    dev = 2.0 / N                            # the truth: exponent -1
    with watch() as w:
        p = scaling(N, dev, "deviation", expect=-0.5)
    assert abs(p + 1.0) < 1e-9
    assert "measured exponent is not the claimed one" in kinds(w)
    assert "-1.0000" in w.findings[0].detail


def test_scaling_confirms_a_right_claim_quietly():
    N = np.array([10.0, 20, 40, 80, 160, 320])
    dev = 2.0 / np.sqrt(N)
    with watch() as w:
        p = scaling(N, dev, "deviation", expect=-0.5)
    assert abs(p + 0.5) < 1e-9
    assert w.clean


def test_scaling_refuses_when_it_is_not_a_power_law():
    N = np.array([10.0, 20, 40, 80, 160, 320])
    y = np.array([1.0, 3.0, 0.4, 9.0, 0.2, 5.0])
    with watch() as w:
        scaling(N, y, "junk")
    assert "not a clean power law" in kinds(w)


def test_scaling_says_so_when_there_is_no_data():
    with watch() as w:
        p = scaling([1.0, 2.0], [1.0, -1.0], "two points")
    assert math.isnan(p)
    assert "cannot fit a slope" in kinds(w)


# -- reporting --------------------------------------------------------------

def test_report_is_clean_and_names_its_own_blind_spot():
    lines = []
    with watch() as w:
        pass
    assert w.report(out=lines.append) is True
    assert "clean" in lines[0]
    assert "typo" in lines[0]                # the report never oversells itself


def test_report_lists_findings_and_counts_repeats():
    lines = []
    with watch() as w:
        for _ in range(3):
            x = np.float64(1e-320) * np.float64(1e-10)
        assert x == 0.0
    assert w.report(out=lines.append) is False
    assert any("x3" in ln for ln in lines)


def test_nesting_restores_numpy_state():
    before = np.geterr()
    with watch():
        with watch():
            pass
    assert np.geterr() == before


# -- THE IMPORTANT ONE ------------------------------------------------------

def test_a_formula_typo_is_invisible():
    """`s` where I meant `s+r`.  Well conditioned, no underflow, no division
    by zero, a perfectly smooth curve, a root that really is a root -- and the
    answer wrong by four orders of magnitude.

    Nothing in this module can see it.  The only thing that saw it in fire 295
    was a second instrument built on a different principle (an exact Monte
    Carlo with no time step).  A clean `dull` report is not a correct result.
    """
    r = 0.7

    def correct(s):
        return 1.0 / np.sqrt(np.float64(s) + r)

    def typo(s):
        return 1.0 / np.sqrt(np.float64(s))          # forgot the + r

    s = np.linspace(1e-4, 5.0, 200)
    with watch() as w:
        a = sane(np.array([correct(v) for v in s]), "correct")
        b = sane(np.array([typo(v) for v in s]), "typo")
        safe_div(np.ones_like(s), s + r, "well conditioned")
    assert w.clean                             # <-- the blind spot, asserted
    assert abs(b[0] / a[0]) > 60               # and the answer is badly wrong
