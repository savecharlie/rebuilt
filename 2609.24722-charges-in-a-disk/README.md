# arXiv:2609.24722 — one hundred thousand charges in a disk

**Lavrov & Nikonov, "Lowest-known Energy Configuration of N = 100 000 Coulomb
Charges in a Disk: Breaking the 10⁵ Barrier"** (21 Sep 2026).

Thirty-one hours on a 24-core workstation produced one number:

```
E_min(100000) = 7.80466624157e9
```

and the paper's headline is that this number sits 2.525×10⁻⁷ away from an
asymptotic expansion whose coefficients were fitted, by other people, to data
that stopped at N = 5000.

**The first two coefficients of that expansion are not fitting parameters.**
Both were available in closed form before any of this was computed — one from
the electrostatics of a charged conducting disk, one from the Madelung energy
of a two-dimensional Wigner crystal. Put together they give

```
E(N) ≈ (π/4) N² − 1.5642653 N^{3/2}
```

with nothing fitted to anything, and that two-term formula reproduces the
31-hour result to **1.9×10⁻⁵**.

This rebuild computes both coefficients from scratch, checks the second
against the fitted value, and asks where the remaining 0.09 % discrepancy
lives.

---

## What is checked, and against what

| file | what it computes | control |
|---|---|---|
| `continuum.py` | k₁ for the uniform and for the equilibrium measure | potential of the conducting disk must be constant; uniform-disk energy must be 8/(3π) |
| `madelung.py` | C_M, Madelung constant of the 2-D triangular OCP | Ewald answer must not depend on the splitting parameter α; u must scale as √n |
| `minimise.py` | my own minimum-energy configurations | N = 2,3,4 closed forms; N = 60,61,92,99 published global minima (arXiv:2609.20777); BLAS kernel vs direct sum; gradient vs finite differences |
| `k2.py` | k₂ = −C_M ∫ρ^{3/2} dA | compared against the fitted k₂ and against the N = 10⁵ datum |
| `fitbasis.py` | where the local-density argument stops being valid | the rim fraction must decay as N^{−1/6} |
| `ladder.py`, `figure.py` | E(N) for 100 ≤ N ≤ 2000, and the running coefficient | |

---

## 1. k₁ = π/4 does not come from a uniform density

The paper says:

> the leading term k₁ = π/4 arises from the average Coulomb repulsion in a
> uniform-density approximation

The uniform-density approximation gives a different number. For unit total
charge spread uniformly over the unit disk,

```
E[uniform]      = 8/(3π) = 0.848826363157      (quadrature: 0.848826362683, rel 5.6e-10)
E[equilibrium]  = π/4    = 0.785398163397      (quadrature: 0.785398086489, rel 9.8e-8)
```

π/4 is the **equilibrium** value — the self-energy of the arcsine measure
ρ(r) = 1/(2π√(1−r²)), the charge distribution on a conducting disk, whose
capacity is C = 2R/π (Thomson, 1867). `continuum.py` confirms it the way a
physicist would rather than by quoting the capacity: it evaluates the
potential of that measure at seven radii and gets π/2 at every one of them, to
ten digits, which is what "conductor" means.

The uniform value is 8.08 % higher. At N = 10⁵ that is a difference of
6.3×10⁸ in the energy — against the 2.0×10³ the paper is celebrating.
The number π/4 is right; the sentence attached to it is not.

## 2. k₂ is the Madelung constant of the local crystal

Where the charge density is n, a triangular lattice sitting in its own
neutralising background has energy per particle −C_M √n. Summing that over the
disk with n = Nρ gives

```
k₂ = −C_M ∫ ρ^{3/2} dA
```

Both factors are computable exactly.

**C_M by Ewald summation** (`madelung.py`), with the splitting parameter as the
control — the split is arbitrary, so a correct implementation cannot depend on
it:

```
alpha = 1.20  u = -1.960515789320
alpha = 1.60  u = -1.960515789320
alpha = 2.00  u = -1.960515789320
alpha = 2.40  u = -1.960515789320
alpha = 2.80  u = -1.960515789320
alpha = 3.20  u = -1.960515789320
spread over alpha = 2.2e-16
```

and u/√n is constant to ten digits over n spanning a factor of 36. The square
lattice gives −1.950132460, i.e. higher, as it must be. Bonsall & Maradudin
(1977) give −1.960516.

**A third road to the same constant, with no physics in it.** `zeta_road.py`
gets C_M out of the analytic continuation of the hexagonal lattice's Epstein
zeta function, which factors as 6 ζ(s/2) L₋₃(s/2) through the Dirichlet
L-function of the quadratic character mod 3:

```
zeta(1/2)             = -1.46035450880959
L_(-3)(1/2)           =  0.480867557696829
zeta_hex(1), covol 1  = -3.92103157863978
u = zeta_hex(1)/2     = -1.96051578931989
Ewald summation       = -1.960515789319892
relative difference   =  1.1e-16
```

Neither side converges anywhere near the point being used — one is a divergent
lattice sum split with a Gaussian, the other is a Dirichlet series continued
past its half-plane — and they agree to the last bit of a double. Controls:
L₋₃ against a direct paired sum (1.7e-17 at s = 4), L₋₃(1+ε) → π/(3√3), and
the lattice sum against its own continuation at s = 6 (5.0e-13).

**The shape functional**, in closed form and by quadrature:

```
∫ ρ_eq^{3/2} dA = 2/√(2π) = 0.797884560803   (quadrature 0.797884560805, rel 3.2e-12)
```

So

```
k₂ = −1.9605157893 × 0.7978845608 = −1.5642652795
```

against the fitted **−1.5628** of Amore & Zárate — a relative difference of
9.4×10⁻⁴.

## 3. What the two free terms predict

```
(π/4)N² − 1.5642653 N^{3/2}   at N = 10⁵   =  7.804515e9
Lavrov & Nikonov, 31 CPU-hours             =  7.804666e9
                                     relative difference 1.9e-5
```

Two constants, one from 1867 and one from 1977, land within 19 parts per
million of a number that took 31 hours to find.

## 4. Where the remaining 0.09 % lives

The local-crystal argument needs the spacing between charges to be small
compared with the distance over which the density changes. At depth d below
the rim, ρ ~ 1/(2π√(2d)), so the spacing is ~N^{−1/2}d^{1/4} and the condition
is d ≫ N^{−2/3}. The argument therefore fails inside an annulus of width
N^{−2/3} — precisely the layer holding the N_b ≈ 2.843 N^{2/3} rim charges
that Amore & Zárate found empirically.

The N^{3/2} weight carried by that annulus is N^{3/2}·(N^{−2/3})^{1/4} =
**N^{4/3}**. `fitbasis.py` measures the fraction and it decays as N^{−1/6}:

```
      N    rim width   tail/total
  1e+02    4.642e-02    5.488e-01
  1e+03    1.000e-02    3.756e-01
  1e+04    2.154e-03    2.561e-01
  1e+05    4.642e-04    1.745e-01
  1e+06    1.000e-04    1.189e-01
  1e+08    4.642e-06    5.520e-02

  ratio per decade: 0.6844, 0.6820, 0.6814, 0.6813      10^(-1/6) = 0.6813
```

Two consequences. First, k₂ = −C_M∫ρ^{3/2} is exact as the coefficient of
N^{3/2} — the failure is confined to a region whose whole contribution is a
lower order. Second, the next term in the expansion is N^{4/3}, and the fitted
basis {N², N^{3/2}, N, N^{1/2}, 1} has no slot for it, so its weight has to go
somewhere. Over 100 ≤ N ≤ 5000 the rim annulus still carries between 55 % and
28 % of the N^{3/2} weight, which is why the fitted k₂ can be pulled 0.09 %
off its limit and still describe the data well.

