# tools

## `dull.py` — the instrument that does the looking

[The ruler is usually the broken thing](../essays/the-ruler-is-usually-the-broken-thing.md)
argues that my instruments have been wrong about twice as often as the papers,
and that the faults are *dull* — a typo, an underflow, a bracket endpoint, a
divide by almost-zero, an inverted guard. They are dull precisely because
there is no reward in looking for them, which is why they survive.

This is the looking, made cheap enough that there is no excuse:

```python
from dull import watch, sane, safe_div, root_ok, scaling

with watch() as w:                                  # numpy FP errors, with the
    ...                                             # CALLER's line number
    sane(profile, "rho(r)")                         # empty / all-nan / constant
    q = safe_div(num, den, "dP/dx")                 # denominator ~ 0
    root_ok(g, root, (a, b), "t_c")                 # root pinned to a dead edge
    p = scaling(Ns, devs, "deviation", expect=-0.5) # MEASURE the exponent
w.report()
```

`watch()` turns on **underflow** reporting, which nobody enables and which
produced a spurious threshold in the rank-resetting rebuild when `Q0**(N-k)`
quietly reached 0.0 at N ≥ 1280. `scaling()` exists because in that same
rebuild I headed a column `dev*√N`, wrote "decays like 1/√N" underneath it,
and did not notice that every row was half the one above — which is 1/N, and
was sitting in the digits. I read my own label instead of the data.

### What it cannot do, asserted rather than hoped

`test_dull.py::test_a_formula_typo_is_invisible` and control 5 of
`dull_validation.py` both build a formula typo — `s` where `s+r` was meant —
that is **wrong by a factor of 84** and passes every check in this module
green. Nothing numerical sees a well-conditioned wrong formula. Only a second
instrument built on a different principle does.

**A clean `dull` report is not a correct result.** It is evidence about the
dull faults only.

### Validation

```
python3 dull_validation.py        # 11/11 controls, against known answers
python3 -m pytest test_dull.py    # 18/18
```

The controls run against the real rebuild code in
[`../2609.23542-rank-resetting/`](../2609.23542-rank-resetting/) with the
historical bug put back: the naive bracket returns `t_c = 0.001` against the
true `0.743904342807` and is caught; the scanned bracket reproduces
`0.743904342807` to 2.2e-13 and is silent.
