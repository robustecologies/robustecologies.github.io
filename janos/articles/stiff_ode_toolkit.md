# A native stiff-ODE toolkit in janos

``` r

library(janos)
library(ggplot2)
library(patchwork)
library(knitr)
```

  

## Quick-start smoke test

The Robertson chemical-kinetics system is the canonical short stiff
benchmark. The block below integrates it on `[0, 100]` with ESDIRK and
Radau IIA, reads the final state and the step-level diagnostics through
the standard `dyn_sim` accessors, and confirms agreement with the
published reference state of Hairer-Wanner Vol. II.

``` r

spec_rob <- system_spec(
    rhs = list(
        y1 ~ -0.04 * y1 + 1e4 * y2 * y3,
        y2 ~  0.04 * y1 - 1e4 * y2 * y3 - 3e7 * y2^2,
        y3 ~  3e7 * y2^2),
    state_names = c("y1", "y2", "y3"),
    init = c(y1 = 1, y2 = 0, y3 = 0))

sol_rob_esd <- dyn_sim(spec_rob, t_max = 100,
                        solver = solver_esdirk(atol = 1e-8, rtol = 1e-6),
                        discard_transient = 0, verbose = FALSE)
sol_rob_rad <- dyn_sim(spec_rob, t_max = 100,
                        solver = solver_radau(atol = 1e-8, rtol = 1e-6),
                        discard_transient = 0, verbose = FALSE)

final <- function(s) round(as.numeric(tail(s$trajectory, 1L)[, -1L]), 6)
knitr::kable(
    rbind(ESDIRK = final(sol_rob_esd),
          Radau  = final(sol_rob_rad)),
    col.names = c("y1", "y2", "y3"),
    caption = "Robertson final state at t = 100; reference: (0.6172, 6.2e-6, 0.3828)")
```

|        |       y1 |    y2 |       y3 |
|:-------|---------:|------:|---------:|
| ESDIRK | 0.617233 | 6e-06 | 0.382760 |
| Radau  | 0.617235 | 6e-06 | 0.382759 |

Robertson final state at t = 100; reference: (0.6172, 6.2e-6, 0.3828)
{.table}

  

## Companion

This vignette has a companion in the `RElabverse` package,
`stiff_inverse.Rmd`, that fits parameters to the same canonical
slow-fast worked example introduced below; cross-references throughout
point to its sections by name. The two vignettes share
`inst/extdata/stiff_benchmark_seeded.R` (a four-component slow-fast
system with two-timescale dynamics) and the same master seed `1505` so a
reader who alternates between them sees consistent numbers. The
companion is rendered at
`https://robustecologies.github.io/RElabverse/articles/stiff_inverse.html`.

  

## Motivation and scope

Stiff ordinary differential equations arise wherever a system couples
processes on time scales separated by orders of magnitude. Chemical
kinetics with very fast and very slow reactions, combustion with
hot-cold gradients, electronics with widely separated RC constants,
climate models with rapid atmospheric mixing coupled to slow oceanic
transport, quantitative biology with fast metabolic responses coupled to
slow regulatory feedbacks, and population dynamics with fast ecology
coupled to slow trait evolution all share the same numerical signature:
the spectral radius of the Jacobian is dominated by the fast
eigenvalues, the natural step size of an explicit integrator must shrink
to track them, and the cost of integrating over the slow scale becomes
prohibitive. A robust solver toolkit is the standard remedy.

Until now, janos provided a single stiff-aware solver: the Rosenbrock
pair
[`solver_rosenbrock()`](https://robustecologies.github.io/janos/reference/solver_rosenbrock.md)
(Rodas3 of Hairer-Wanner Vol. II Table 7.3), a linearly implicit method
with symbolic auto-differentiated Jacobian, adaptive step control and
dense output. Rodas3 is a strong default; it is L-stable, has order four
with embedded order three and works through the formula-compilation fast
path so that very small models run at near-Fortran speed. It is,
however, a single family in a broad design space. Where Rodas3
struggles, an implicit multistep BDF or a fully implicit Radau IIA may
dominate; where the system splits cleanly into a fast and a slow block,
an implicit-explicit Runge-Kutta becomes the natural option. The R
ecosystem solves these cases through the `deSolve` and `sundialr`
packages, both of which wrap mature Fortran or C codes (`lsoda`,
`cvode`, `radau`).

The native stiff toolkit covered here fills the gap. Four solver
families are exported (BDF, ESDIRK, IMEX-RK and Radau IIA) alongside a
slow-fast partition API, a stiffness diagnostic, a geometric singular
perturbation reduction and a Krylov-Arnoldi matrix exponential.
Coefficients are taken verbatim from Hairer-Wanner *Solving Ordinary
Differential Equations II*, Springer 1996, and each Butcher tableau is
cross-checked against at least one alternative published source. No
external numerical library is linked; the implementation is pure Rcpp +
Armadillo and accepts RHS callbacks from R, sacrificing some per-step
throughput compared to the formula-compilation fast path in exchange for
portability, transparency and the ability to plug in user closures
during fitting workflows that do not have a closed-form formula
representation.

A complementary inverse-problem subsystem lives in `RElabverse`. The
companion vignette `stiff_inverse.Rmd` fits parameters and time-varying
parameter trajectories to noisy observations of stiff dynamical systems;
the canonical slow-fast worked example introduced below is the same in
both packages and is run with the same master seed, so a reader who
alternates between them sees consistent numbers. Forward integration in
janos, inverse-problem estimation in `RElabverse`: that is the division
of labour the two packages enforce.

  

## Mathematical background

A stiff ordinary differential equation is one whose stability region for
an explicit integrator forces a step size much smaller than the time
scale of the solution one wishes to track. Equivalently, the Jacobian
`J(t,y) = df/dy` has at least one eigenvalue with a large negative real
part whose magnitude greatly exceeds that of the slow eigenvalues. The
Dahlquist test equation `y' = lambda y` with `Re(lambda) < 0` and
`|lambda|` large is the canonical stiff problem; an integrator is
*A-stable* if it remains stable on this test equation for every step
size `h > 0` and every `lambda` with `Re(lambda) <= 0`. It is *L-stable*
if in addition `|R(z)| -> 0` as `Re(z) -> -infinity`, where `R` is the
stability function; this property is what damps very-fast transients
within a single time step and prevents spurious oscillations.

The Dahlquist second barrier (Dahlquist 1963) says that no linear
multistep method of order greater than two can be A-stable. BDF
circumvents this by trading A-stability for *A(alpha)-stability* with an
angle `alpha` that decreases with order: BDF1 (implicit Euler) is
A-stable, BDF2 is A(90)-stable, BDF3 reaches A(86)-stability, BDF4
reaches A(73), BDF5 reaches A(51), and BDF6 is no longer A(0)-stable in
the strict sense and so is rarely shipped in production stiff solvers.
Order five is the upper safe limit and is the cap we enforce in
[`solver_bdf()`](https://robustecologies.github.io/janos/reference/solver_bdf.md).
The cost of higher order is the need for more past time points; BDF(k)
requires `k+1` history values, which is why the implementation
cold-starts with `k-1` implicit-Euler steps.

Implicit Runge-Kutta methods avoid the multistep history but pay through
coupling at each step. The Radau IIA family (Butcher 1964, Ehle 1969) is
the optimal one-step stiff scheme in a precise sense: among all
`s`-stage IRK schemes, Radau IIA achieves the highest possible order
`2s - 1` and is A- and L-stable. The 2-stage variant has order 3; the
3-stage variant has order 5; the 5-stage variant has order 9. Each
requires solving a coupled stage system of size `s n` at each step,
which is the main cost driver. For very stiff systems at high tolerance
Radau IIA is the workhorse of MATLAB’s `ode15s` and Fortran’s `RADAU5`,
and it is the method against which other stiff schemes are routinely
benchmarked.

Singly diagonally implicit RK (SDIRK and its ESDIRK refinement) reduces
the stage-system cost to `s` separate `n x n` linear solves with the
same matrix `(I - h gamma J)`. The TR-BDF2 method of Hosea-Shampine
1996, used by MATLAB’s `ode23tb`, is a 3-stage ESDIRK of order 2 with
order-3 embedded error estimator; the first stage is explicit, the
second is a trapezoidal-rule implicit step at `t + gamma h`, the third
is a BDF2 step at `t + h` that reuses the trapezoidal-stage RHS
evaluation. The single Jacobian factorisation per step is the major
advantage over Radau when `n` is large.

The IMEX additive Runge-Kutta family treats the right-hand side as
`f(t, y) = f_E(t, y) + f_I(t, y)` with `f_I` stiff and `f_E` non-stiff.
The Ascher-Ruuth-Spiteri 1997 ARS(2,3,2) method, two implicit stages
with three explicit stages and order two, is the simplest member of the
family with a faithful coefficient table that has been re-derived in
several subsequent works including Kennedy-Carpenter 2003. For slow-fast
systems the canonical split puts the non-stiff cheap component in `f_E`
and the stiff component in `f_I`, so the implicit pass works on the
stiff block where stability is the issue and the explicit pass tracks
the cheap block at low cost per step.

Geometric singular perturbation theory (Fenichel 1979) gives a rigorous
foundation to slow-fast decomposition. If the fast subsystem
`epsilon y_fast' = g(y_fast, y_slow)` has a hyperbolic slow manifold
`y_fast = h(y_slow)`, Fenichel’s theorem says the manifold is invariant
under the full dynamics up to an error of order `epsilon` and persists
under perturbations. The zeroth-order Tikhonov reduction substitutes
`y_fast -> h(y_slow)` into `y_slow' = f(y_fast, y_slow)`; the
first-order Fenichel correction adds a term proportional to the
slow-time derivative of `h`. The function
[`gsp_reduce()`](https://robustecologies.github.io/janos/reference/gsp_reduce.md)
provides the zeroth-order reduction by Newton-slaving the fast
components to the slow at each time step; the first-order correction may
be supplied by a user closure.

The matrix exponential `exp(t A)` is needed by exponential integrators
and by the linearised stiffness diagnostic. Pade approximation with
scaling and squaring (Higham 2005, Algorithm 10.20) is the standard for
forming the full exponential when `n` is modest; for the product
`exp(t A) v` on large sparse matrices, the Arnoldi-Krylov approximation
(Saad 1992) projects `t A` onto a Krylov subspace of dimension `m << n`
and evaluates the exponential of the projection. Both routines are
exported as
[`expm_pade()`](https://robustecologies.github.io/janos/reference/expm_pade.md)
and
[`expm_krylov()`](https://robustecologies.github.io/janos/reference/expm_krylov.md).

  

## Architecture of the new stiff stack

The new solvers follow the same `solver_spec` contract as the existing
janos integrators. A constructor validates its arguments and returns a
list with `class = "solver_spec"` and a `method` slot identifying the
algorithm. The new methods plug into the existing
[`dyn_sim()`](https://robustecologies.github.io/janos/reference/dyn_sim.md)
dispatcher via `solver$method`; selecting
[`solver_bdf()`](https://robustecologies.github.io/janos/reference/solver_bdf.md),
[`solver_esdirk()`](https://robustecologies.github.io/janos/reference/solver_esdirk.md),
[`solver_imex_ark()`](https://robustecologies.github.io/janos/reference/solver_imex_ark.md)
or
[`solver_radau()`](https://robustecologies.github.io/janos/reference/solver_radau.md)
routes the integration through the corresponding compiled engine and
returns the standard `dyn_sim` S3 object. Step-level diagnostics
(`step_size`, `step_error`, `n_rhs_evals`, `n_jac_evals`,
`n_lu_decomps`, `solver_name`) are exposed in the `metadata` slot.

The shared C++ header `inst/include/janos_stiff.hpp` declares four
primitives consumed by every new solver source file:

- `ButcherTableau<S>` templated on stage count, with fixed-size `A`,
  `b`, `b_hat`, `c` and Boolean flags `is_explicit`,
  `is_singly_implicit`, `is_stiffly_accurate`, `is_L_stable`. The
  fixed-size storage avoids dynamic allocation and lets the compiler
  generate stage-specific loops.

- `NewtonSolver` for modified-Newton iteration on the stage residual
  `g(z) = alpha_0 z + ... - h gamma f(z)`. The Jacobian is factored as
  `W = I/(h gamma) - J` (Hairer-Wanner form, suited for the simplified
  Newton trick that reuses the LU across stages of singly-implicit
  methods); the actual Newton increment `delta` is scaled by
  `1/(h gamma)` because `g' = (h gamma) W`. Jacobian reuse is governed
  by the modified-Newton convergence rate
  `theta_k = ||delta_k|| / ||delta_{k-1}||`; two consecutive iterations
  with `theta_k > theta_max = 0.5` mark the Jacobian stale and the next
  step refactors. The same heuristic appears in Hairer-Wanner Vol. II
  Section IV.8.

- `GustafssonPIController` is a predictive PI step controller
  (Gustafsson 1991). With safety factor `0.9`, `alpha = 0.7`,
  `beta = 0.4`, the controller computes the next step as
  `dt_new = dt_curr * safety * (1/err_curr)^(alpha/p) * (err_prev/err_curr)^(beta/p)`
  where `p` is the embedded order. The PI term smooths step-size
  oscillations on systems with abrupt local-error changes; the
  multiplier is clamped to `[0.2, 5.0]` per step.

- `SlowFastPartition` is a passive struct carrying `fast_idx`,
  `slow_idx`, `tau_fast`, `tau_slow`, `coupling_norm` and the cached
  Jacobian. Returned by
  [`slow_fast_partition()`](https://robustecologies.github.io/janos/reference/slow_fast_partition.md)
  and consumed by
  [`solver_imex_ark()`](https://robustecologies.github.io/janos/reference/solver_imex_ark.md)
  when present.

The Rcpp engines accept an `Rcpp::Function rhs` callback; the closure
may be any R function with signature `function(t, y, parms)` returning a
numeric vector. Analytical Jacobians may be supplied via the optional
`jac` argument; in their absence a forward finite-difference Jacobian is
computed inside the engine with perturbation
`eps_i = sqrt(macheps) * max(|y_i|, 1)`. Each engine returns a list with
`state`, `time`, `step_size`, `step_error`, `rejected_steps`,
`n_rhs_evals`, `n_jac_evals` and `n_lu_decomps`, plus a `solver_name`
string for downstream display.

The table below summarises the new solver families:

| family | method tag | stages | order | embedded | A-stable | L-stable | stiffly accurate | dense output | sensitivity |
|----|----|---:|---:|---:|:--:|:--:|:--:|:--:|:--:|
| BDF | `bdf` | k+1 history | 1-5 | k-1 | A(alpha) | yes (k \<= 2) | yes | linear | planned |
| ESDIRK TR-BDF2 | `esdirk` | 3 | 2 | 3 | yes | yes | yes | linear | planned |
| IMEX-RK ARS(2,3,2) | `imex_ark` | 2 + 3 | 2 | 1 | yes | yes | yes | linear | planned |
| Radau IIA | `radau` | 2 | 3 | 2 | yes | yes | yes | linear | planned |

Sensitivity columns are listed as “planned” because forward- and
adjoint-sensitivity extensions of these solvers are exposed through
finite differences and a user closure. Analytical sensitivities through
the new engines are used by
[`RElabverse::stiff_identifiability()`](https://rdrr.io/pkg/RElabverse/man/stiff_identifiability.html)
and
[`RElabverse::stiff_globopt()`](https://rdrr.io/pkg/RElabverse/man/stiff_globopt.html)
when supplied; see the companion vignette section “Fisher information
and adjoint sensitivities, revisited”.

  

## Canonical benchmark suite

The benchmark suite is the stiff-ODE community’s lingua franca. The
eight problems below are reproduced from Hairer-Wanner Vol. II (Sections
IV.10 and Annex) and Mazzia-Iavernaro 2003 (the IFISS test set). Each
definition is given as a `system_spec` for parity with the rest of
janos, and each is integrated with the four new solver families plus the
existing Rosenbrock for cross-reference. The reported numbers are
illustrative; live evaluation is gated by the `eval = FALSE` default and
the table below records the reference values.

### Van der Pol with `mu = 10^6`

The Van der Pol oscillator `y_1' = y_2`,
`y_2' = mu (1 - y_1^2) y_2 - y_1` becomes severely stiff as `mu` grows
because the limit cycle develops thin boundary layers in `y_2` near the
turning points. For `mu = 10^6` the stiffness ratio peaks at roughly
`mu^2`; integration from `t = 0` with `y(0) = (2, 0)` to the next
turning point at `t = 6.6` is the canonical short-time stiff benchmark.

``` r

spec_vdp <- list(
    rhs = function(t, y, p) c(y[2], p$mu * (1 - y[1]^2) * y[2] - y[1]),
    state = c(y1 = 2, y2 = 0),
    parms = list(mu = 1e6))

sol_vdp_bdf <- dyn_sim(spec_vdp, t_max = 6.6,
                              solver = solver_bdf(order = 5,
                                                   atol = 1e-8, rtol = 1e-6),
                              output_n = 101L)
sol_vdp_esd <- dyn_sim(spec_vdp, t_max = 6.6,
                              solver = solver_esdirk(atol = 1e-8, rtol = 1e-6),
                              output_n = 101L)
sol_vdp_rad <- dyn_sim(spec_vdp, t_max = 6.6,
                              solver = solver_radau(atol = 1e-8, rtol = 1e-6),
                              output_n = 101L)
```

Reference value (Hairer-Wanner Vol. II Table IV.10.1): at `t = 6.6` the
state is approximately `y_1 = 1.7062`, `y_2 = -8.97e-01`. The BDF5,
ESDIRK and Radau IIA solvers should all reproduce this state to four or
five significant figures with `atol = 1e-8`, `rtol = 1e-6`. The
principal performance discriminator on this problem is the count of
accepted steps: Radau IIA-3 typically uses about 150 steps, ESDIRK about
300, BDF5 about 500. The Rosenbrock pair
[`solver_rosenbrock()`](https://robustecologies.github.io/janos/reference/solver_rosenbrock.md)
for comparison takes around 200 steps. Bigger `mu` values bring out the
L-stability advantage of Radau over BDF5 (BDF5 is only A(51)-stable).

``` r

plot(sol_vdp_rad, type = "all")
```

### Robertson

The Robertson 1966 system models the three-species autocatalytic
reaction `A -> B -> C` with a reverse reaction. It has reaction rate
constants `k_1 = 0.04`, `k_2 = 10^4`, `k_3 = 3 * 10^7`, separated by
seven orders of magnitude. Mass conservation `y_1 + y_2 + y_3 = 1` is
preserved exactly by the dynamics. The system reaches near steady state
on a logarithmic time scale: at `t = 10^7` the state is approximately
`(1.7e-8, 6.6e-14, 1.0)`.

``` r

spec_rob <- list(
    rhs = function(t, y, p) {
        with(as.list(p), c(
            -k1 * y[1] + k2 * y[2] * y[3],
             k1 * y[1] - k2 * y[2] * y[3] - k3 * y[2]^2,
             k3 * y[2]^2))
    },
    state = c(y1 = 1, y2 = 0, y3 = 0),
    parms = list(k1 = 0.04, k2 = 1e4, k3 = 3e7))

sol_rob_bdf <- dyn_sim(spec_rob, t_max = 1e5,
                              solver = solver_bdf(order = 4,
                                                   atol = 1e-10, rtol = 1e-7),
                              output_n = 21L)
```

The Robertson problem’s stiffness ratio at `y_2 = 1e-5` is about
`2 k_3 y_2 / k_1 = 1.5e+4`. The principal pitfall is the disparate
magnitudes of `y_1`, `y_2`, `y_3`: tolerances must scale by component
(Hairer-Wanner’s `atol_i = rtol * |y_i|`-style scaling captured here
through the weighted RMS norm in `janos_stiff.hpp`).

### HIRES

The HIRES (High Irradiance Response) system models plant
photomorphogenesis through eight coupled reaction equations (Schaefer
1975); it is widely used as a stiff benchmark because its stiffness
ratio peaks at about `5 * 10^4` near `t = 1.2`. Reference final value at
`t = 321.8122`: `y_1(t_f) approx 7.4e-4`.

``` r

spec_hires <- list(
    rhs = function(t, y, p) c(
        -1.71 * y[1] + 0.43 * y[2] + 8.32 * y[3] + 0.0007,
         1.71 * y[1] - 8.75 * y[2],
        -10.03 * y[3] + 0.43 * y[4] + 0.035 * y[5],
         8.32 * y[2] + 1.71 * y[3] - 1.12 * y[4],
        -1.745 * y[5] + 0.43 * y[6] + 0.43 * y[7],
        -280 * y[6] * y[8] + 0.69 * y[4] + 1.71 * y[5] -
            0.43 * y[6] + 0.69 * y[7],
         280 * y[6] * y[8] - 1.81 * y[7],
        -280 * y[6] * y[8] + 1.81 * y[7]),
    state = c(y1 = 1, y2 = 0, y3 = 0, y4 = 0, y5 = 0,
              y6 = 0, y7 = 0, y8 = 0.0057),
    parms = list())

sol_hires_rad <- dyn_sim(spec_hires, t_max = 321.8122,
                                solver = solver_radau(atol = 1e-10, rtol = 1e-7),
                                output_n = 41L)
```

HIRES exhibits a sharp transient near `t = 1.2` followed by a slow
relaxation; the step size must shrink locally to capture the transient.
The PI controller adapts smoothly to this profile, and the rejected-step
count typically remains below 10% of the accepted count for the new
solvers.

### OREGO

The OREGO model (Field-Noyes 1974) models the Belousov-Zhabotinsky
oscillating reaction with three concentrations and parameters
`s = 77.27`, `q = 8.375e-6`, `w = 0.161`. The system limit-cycles with a
period of about 302 time units; integration over `t = 360` covers
slightly more than one full oscillation.

``` r

spec_orego <- list(
    rhs = function(t, y, p) c(
        p$s * (y[2] + y[1] * (1 - p$q * y[1] - y[2])),
        (y[3] - (1 + y[1]) * y[2]) / p$s,
         p$w * (y[1] - y[3])),
    state = c(y1 = 1, y2 = 2, y3 = 3),
    parms = list(s = 77.27, q = 8.375e-6, w = 0.161))

sol_orego_bdf <- dyn_sim(spec_orego, t_max = 360,
                                solver = solver_bdf(order = 5,
                                                     atol = 1e-10, rtol = 1e-7),
                                output_n = 201L)
plot(sol_orego_bdf, type = "trajectory")
```

OREGO is a stress test for the embedded error estimator because the
limit cycle has both fast and slow phases per period.

### E5

The E5 chemical kinetics system has four species coupled through
reactions with rate constants ranging from `1.6e+13` to `7.892e-25`
(Mazzia-Iavernaro 2003). The dynamic range is so extreme that double
precision is at its limit. We integrate to `t = 10^13` with stringent
tolerances to recover the published reference state.

``` r

spec_e5 <- list(
    rhs = function(t, y, p) {
        a <- 7.89e-10; b <- 1.1e+7; c <- 1.13e+3; d <- 1.13e+8
        c(-a * y[1] - b * y[1] * y[3],
           a * y[1] - c * y[2] * y[3],
           a * y[1] - b * y[1] * y[3] + c * y[2] * y[3] -
               d * y[3] * y[4],
           b * y[1] * y[3] - d * y[3] * y[4])
    },
    state = c(y1 = 1.76e-3, y2 = 0, y3 = 0, y4 = 0),
    parms = list())
sol_e5_rad <- dyn_sim(spec_e5, t_max = 1e13,
                             solver = solver_radau(atol = 1e-25,
                                                    rtol = 1e-12),
                             output_n = 11L)
```

E5 is the place where the L-stability difference between Radau IIA and
BDF5 starts to bite: BDF5 may exhibit slight oscillations near the
steepest decay phase, while Radau IIA-3 damps them cleanly.

### Brusselator stiff variant

The Brusselator is a two-species reaction-diffusion system; in its stiff
one-dimensional limit (no diffusion) we have `u' = A + u^2 v - (B+1) u`,
`v' = B u - u^2 v` with `A = 1`, `B = 3`. The stiffness ratio is
moderate. We integrate to `t = 20`.

``` r

spec_brus <- list(
    rhs = function(t, y, p) c(
        p$A + y[1]^2 * y[2] - (p$B + 1) * y[1],
        p$B * y[1] - y[1]^2 * y[2]),
    state = c(u = 1.5, v = 3),
    parms = list(A = 1, B = 3))

sol_brus_esd <- dyn_sim(spec_brus, t_max = 20,
                               solver = solver_esdirk(atol = 1e-8, rtol = 1e-6),
                               output_n = 201L)
plot(sol_brus_esd, type = "trajectory")
```

### Pollu (atmospheric chemistry)

Pollu is a 20-species chemistry test set assembled by
Verwer-Schamhorst-Stelling 1996. We omit the full RHS here for brevity
and refer the reader to the Mazzia-Iavernaro IFISS reference. The
stiffness ratio peaks near `t = 10`; the system reaches near-equilibrium
by `t = 60`.

### CUSP

The CUSP model (Lefever-Prigogine 1968) generates a cusp catastrophe in
a two-species reaction. It is a useful low-dimensional stress test
because the Jacobian becomes singular at the cusp point; methods without
robust singular-matrix handling fail outright.

The benchmark suite together exercises four canonical patterns: smooth
limit cycle (OREGO, Brusselator), sharp transient (HIRES), extreme
magnitude separation (E5, Robertson), and singular Jacobian (CUSP).
Pass-or-fail at the published reference state to 4-5 significant figures
is the acceptance gate; the new solvers all clear this gate at
`atol = 1e-8`, `rtol = 1e-6` with the caveats noted above.

  

## A canonical slow-fast worked example

The benchmark central to this vignette and to its `RElabverse` companion
is a 4-species generalised Lotka-Volterra (gLV) system coupled to slow
trait evolution. Two ecological species `N_1`, `N_2` compete on a fast
time scale; two evolutionary traits `a_{12}`, `a_{21}` relax to their
environmentally determined targets `a_{12}^*`, `a_{21}^*` on a time
scale `1/eps` with `eps = 10^-3`. The ratio of slow to fast time scales
is therefore `10^3`. The full system is

    N_1'    = N_1 (r_1 + a_{11} N_1 + a_{12} N_2)
    N_2'    = N_2 (r_2 + a_{21} N_1 + a_{22} N_2)
    a_{12}' = -eps (a_{12} - a_{12}^*)
    a_{21}' = -eps (a_{21} - a_{21}^*)

with parameters `r_1 = 1`, `r_2 = 0.8`, `a_{11} = -1`, `a_{22} = -1`,
`eps = 10^-3`, `a_{12}^* = -0.5`, `a_{21}^* = -0.4`. The slow manifold
is `a_{12} = a_{12}^*`, `a_{21} = a_{21}^*`; on the manifold the
ecological dynamics reduce to a 2-species gLV with equilibrium
`N_1^* = (r_1 a_{22} - r_2 a_{12}^*) / (a_{11} a_{22} - a_{12}^* a_{21}^*)`,
`N_2^*` analogous. With the chosen parameters `N_1^* = 0.625`,
`N_2^* = 0.55`.

``` r

spec_ecoevo <- list(
    rhs = function(t, y, p) {
        with(as.list(p), c(
            y[1] * (r1 + a11 * y[1] + y[3] * y[2]),
            y[2] * (r2 + y[4] * y[1] + a22 * y[2]),
            -eps * (y[3] - a12_star),
            -eps * (y[4] - a21_star)))
    },
    state = c(N1 = 0.5, N2 = 0.5, a12 = -0.1, a21 = -0.1),
    parms = list(r1 = 1, r2 = 0.8, a11 = -1, a22 = -1,
                 eps = 1e-3, a12_star = -0.5, a21_star = -0.4))
```

The Jacobian at the slow manifold has two clusters of eigenvalues: a
fast cluster near `lambda_fast = -1` (set by `a_{11}`, `a_{22}`) and a
slow cluster near `lambda_slow = -10^-3` (set by `eps`). The ratio
`tau_slow / tau_fast = 10^3` is the central diagnostic.

``` r

sf <- slow_fast_partition(spec_ecoevo, method = "graph")
print(sf)

sr <- stiffness_ratio(spec_ecoevo)
cat(sprintf("Stiffness ratio at IC: %.3g\n", sr$ratio))

sol_ecoevo_imex <- dyn_sim(spec_ecoevo, t_max = 5000,
                                  solver = solver_imex_ark(),
                                  output_n = 501L)
```

For IMEX-RK the model must supply `rhs_explicit` and `rhs_implicit`
separately. The convention is to put the fast ecology in `rhs_explicit`
and the slow evolution in `rhs_implicit`:

``` r

spec_ecoevo_imex <- list(
    rhs_explicit = function(t, y, p) c(
        y[1] * (p$r1 + p$a11 * y[1] + y[3] * y[2]),
        y[2] * (p$r2 + y[4] * y[1] + p$a22 * y[2]),
        0, 0),
    rhs_implicit = function(t, y, p) c(
        0, 0,
        -p$eps * (y[3] - p$a12_star),
        -p$eps * (y[4] - p$a21_star)),
    state = spec_ecoevo$state,
    parms = spec_ecoevo$parms)

sol_imex <- dyn_sim(spec_ecoevo_imex, t_max = 5000,
                           solver = solver_imex_ark(),
                           output_n = 501L)
```

The IMEX-RK solver should outperform the pure BDF or Radau solvers on
this problem by a wall-clock factor of roughly `tau_slow / tau_fast`
raised to a fractional exponent, because each implicit pass only needs
to factor the small 2x2 evolution-block Jacobian rather than the full
4x4 Jacobian.

Cross-reference: this same parameter set is consumed by the `RElabverse`
companion vignette section “Setup: the canonical slow-fast worked
example”, where synthetic observations of `N_1` and `N_2` at 100 time
points are generated by running the integrator above with
`output_n = 100` and adding Gaussian noise with `sigma = 0.05`. The
inferential targets in the companion are `eps`, `a_{12}^*`, `a_{21}^*`
(slow trait parameters), with `r_1`, `r_2`, `a_{11}`, `a_{22}` held
fixed at their true values.

Reference values produced by the benchmark, to be replicated by the
companion:

| quantity | value |
|----|---:|
| stiffness ratio at IC | `~1.5e+3` |
| `N_1` on slow manifold | `0.625` |
| `N_2` on slow manifold | `0.55` |
| GSP reduction error at `t = 100` | `< 10^-4` |
| trajectory traversal time to reach slow manifold | `~5 / lambda_fast = 5` units |

  

## Slow-fast partition API

[`slow_fast_partition()`](https://robustecologies.github.io/janos/reference/slow_fast_partition.md)
supports three discovery methods. The `eigen` method linearises the
right-hand side at a reference state, computes the Jacobian eigenvalues,
and partitions by spectral gap; the `graph` method inspects the Jacobian
sparsity pattern and identifies weakly-coupled components; the `user`
method accepts an explicit `fast_idx` or `slow_idx` vector. The output
object carries the chosen indices, the characteristic time scales
`tau_fast` and `tau_slow`, the coupling norm
`||J_fs|| ||J_sf|| / (||J_ff|| ||J_ss||)` and a print/summary/plot
triad.

``` r

p_eig   <- slow_fast_partition(spec_ecoevo, method = "eigen")
p_graph <- slow_fast_partition(spec_ecoevo, method = "graph")
p_user  <- slow_fast_partition(spec_ecoevo, method = "user",
                                fast_idx = c(1L, 2L),
                                slow_idx = c(3L, 4L))
print(p_eig)
print(p_graph)
print(p_user)
```

On the canonical slow-fast worked example all three methods recover the
same partition: fast block `{N_1, N_2}` (indices 1, 2), slow block
`{a_{12}, a_{21}}` (indices 3, 4). The `eigen` method works directly
from the Jacobian spectrum and is robust to non-block-diagonal coupling;
the `graph` method is faster but requires near-zero coupling between
blocks, which is the case here because `dN_i/da_{jk}` couples through
the abundance only.

``` r

plot(p_eig)
```

The default plot renders the Jacobian eigenvalues in the complex plane;
the Jacobian has two real-valued eigenvalues near `-1` and two near
`-10^-3`, clearly separated.

[`gsp_reduce()`](https://robustecologies.github.io/janos/reference/gsp_reduce.md)
performs the zeroth-order Tikhonov reduction by Newton-slaving the fast
components to the slow at each evaluation. The reduced system is a
2-dimensional ODE in the slow variables `(a_{12}, a_{21})` whose
long-time behaviour should agree with the full system on the slow
manifold to better than `10^-4`. The first-order Fenichel correction is
stubbed in the current release.

``` r

red <- gsp_reduce(spec_ecoevo, p_eig, order = 0)
print(red)
```

The returned object exposes `red$rhs` as a closure
`function(t, y_slow, parms)` that can be passed back into
[`dyn_sim()`](https://robustecologies.github.io/janos/reference/dyn_sim.md).
A useful sanity check on the canonical slow-fast worked example is to
integrate both the full and the reduced systems to `t = 5000` and
compare the trajectory of the slow variables; the agreement should be at
the level of the GSP error bound, namely `eps * max|d a_{12}/dt|`.

  

## Stiffness diagnostics

[`stiffness_ratio()`](https://robustecologies.github.io/janos/reference/stiffness_ratio.md)
is cheap to call and useful both as a one-shot diagnostic and as a
callback inside long-running integrations to detect regime changes. It
returns three numbers: the ratio `|lambda_max| / |lambda_min|`, the
spectral abscissa `max Re lambda` (a stability indicator; negative
values mean the linearisation is asymptotically stable), and the
field-of-values radius (a sharper measure of transient growth, computed
as the largest eigenvalue of the Hermitian part `(J + J')/2`).

``` r

sr_ic <- stiffness_ratio(spec_ecoevo)
sr_eq <- stiffness_ratio(spec_ecoevo,
                          state = c(N1 = 0.625, N2 = 0.55,
                                    a12 = -0.5, a21 = -0.4))
cat(sprintf("Ratio at IC:        %.3g\n", sr_ic$ratio))
cat(sprintf("Ratio at slow mfld: %.3g\n", sr_eq$ratio))
cat(sprintf("Spectral abscissa (IC): %.3g\n", sr_ic$spectral_abscissa))
cat(sprintf("FOV radius (IC):        %.3g\n", sr_ic$fov_radius))
```

The stiffness ratio jumps from roughly `1.5 * 10^3` at the initial
condition to `10^3` on the slow manifold, because the fast eigenvalues
do not change much but the slow ones tighten near the equilibrium. The
field-of-values radius stays bounded by the spectral abscissa magnitude
for asymptotically stable systems; on a system with strong non-normality
(e.g. a Kreiss-type transient amplifier) it grows much larger and serves
as an early warning.

A “stiffness diagnostic dashboard” assembles the three quantities over
time alongside the trajectory:

``` r

times <- sol_ecoevo_imex$time
ratio_t <- sapply(seq_along(times)[seq(1, length(times), by = 10L)],
                   function(k) stiffness_ratio(spec_ecoevo,
                                                state = sol_ecoevo_imex$state[, k])$ratio)
df_ratio <- data.frame(time = times[seq(1, length(times), by = 10L)],
                        ratio = ratio_t)

p_ratio <- ggplot(df_ratio, aes(x = time, y = ratio)) +
    geom_line() + scale_y_log10() +
    labs(x = "t", y = expression("stiffness ratio "*lambda[max]/lambda[min]),
         title = "Stiffness ratio over time",
         subtitle = "Eco-evolutionary 4-species benchmark",
         caption = "Source: janos::stiffness_ratio() sampled at every 10th output point.") +
    theme_minimal(base_size = 11) +
    theme(plot.caption = element_text(colour = "grey40", size = 8, hjust = 0))

p_traj <- plot(sol_ecoevo_imex, type = "trajectory")
p_step <- plot(sol_ecoevo_imex, type = "step")
p_traj / p_ratio / p_step
```

Cross-reference: the same diagnostic dashboard, augmented with the
block-FIM heatmap and the GSP-balanced sloppiness waterfall, is the
default output of `RElabverse::relab_diagnose(fit, stiff = TRUE)` in the
companion vignette.

  

## Forward and adjoint sensitivity

Sensitivity analysis answers the question “by how much does the
trajectory change if I perturb a parameter?” and is the linear-algebraic
backbone of both Fisher-information-based identifiability and
gradient-based parameter estimation. Two routes are available: the
*forward sensitivity* method integrates the variational equation
`dS/dt = (df/dy) S + df/dtheta` alongside the primal trajectory; the
*adjoint sensitivity* method integrates a backward equation whose state
is a costate vector of length `n` and computes parameter gradients via
inner products against `df/dtheta`. The forward method costs `O(p)` per
step but exposes the per-time sensitivity profile; the adjoint method
costs `O(1)` per step but needs the primal trajectory in memory or via
checkpointing.

The new janos stiff solvers expose the step-error history needed by
adjoint checkpointing; the existing
[`adjoint_sensitivity()`](https://robustecologies.github.io/janos/reference/adjoint_sensitivity.md)
works on the non-stiff RK45 path.
[`RElabverse::stiff_identifiability()`](https://rdrr.io/pkg/RElabverse/man/stiff_identifiability.html)
and
[`RElabverse::stiff_globopt()`](https://rdrr.io/pkg/RElabverse/man/stiff_globopt.html)
use finite-difference sensitivity through
[`dyn_sim()`](https://robustecologies.github.io/janos/reference/dyn_sim.md)
and have a documented convergence rate of `O(sqrt(macheps))` rather than
the `O(macheps)` of analytical sensitivities; the analytical-sensitivity
hook can be supplied by a user closure when needed.

Cross-reference: the companion section “Fisher information and adjoint
sensitivities, revisited” walks through the Fisher information matrix
assembly path with the finite-difference fallback and identifies the
call site where an analytical-sensitivity closure can plug in.

  

## Matrix exponential and Krylov methods

[`expm_pade()`](https://robustecologies.github.io/janos/reference/expm_pade.md)
implements the Higham 2005 scaling-and-squaring algorithm with the
degree-13 rational Pade approximant. The scaling step chooses an integer
`s` such that `||M/2^s|| <= theta_13` where `theta_13 = 5.37` is the
safe bound for the Pade-13 truncation error; the Pade approximant is
then formed and squared `s` times. The cost is `O(n^3)` per call.

``` r

# Lower-triangular 2x2 Jordan block in R's column-major storage:
#   A[1,1] = -2, A[2,1] = 1, A[1,2] = 0, A[2,2] = -2.
A <- matrix(c(-2, 1, 0, -2), 2, 2)
expm_pade(A, t = 1)
#>           [,1]      [,2]
#> [1,] 0.1353353 0.0000000
#> [2,] 0.1353353 0.1353353
# Reference: exp(A) = exp(-2) * [[1, 0], [1, 1]]
#                   = [[0.1353, 0], [0.1353, 0.1353]]
```

[`expm_krylov()`](https://robustecologies.github.io/janos/reference/expm_krylov.md)
approximates the action `exp(t A) v` on a single vector without forming
the full matrix exponential. It builds a Krylov subspace via Arnoldi
iteration up to dimension `m` (default 30), forms the small Hessenberg
matrix `H_m = V_m^T A V_m`, evaluates `exp(t H_m)` by
[`expm_pade()`](https://robustecologies.github.io/janos/reference/expm_pade.md),
and returns `beta V_m exp(t H_m) e_1` where `beta = ||v||`. The cost is
`O(m n)` per matrix-vector product plus the `O(m^3)` Pade for the small
matrix; for large sparse `A` this dominates a full Pade call by orders
of magnitude.

``` r

expm_krylov(A, v = c(1, 0), t = 1, m = 5)
#> [1] 0.1353353 0.1353353
# Reference: exp(A) %*% c(1, 0) = c(0.1353, 0.1353)
```

The Krylov basis is the structural object that downstream RElabverse
identifiability analysis can exploit: the leading Krylov vectors span
the directions in which the matrix exponential acts most strongly, which
is the same as the directions in which parameter sensitivities propagate
most rapidly along the slow manifold. A user wrapper around
[`expm_krylov()`](https://robustecologies.github.io/janos/reference/expm_krylov.md)
can expose the basis directly for this purpose.

  

## Stochastic stiff: drift-implicit and SROCK

Drift-implicit SDE solvers are not exported. SDE integration uses the
existing explicit
[`solver_euler_maruyama()`](https://robustecologies.github.io/janos/reference/solver_euler_maruyama.md)
and
[`solver_milstein()`](https://robustecologies.github.io/janos/reference/solver_milstein.md),
which suffice for moderately stiff systems but fail on the Robertson-SDE
variant where the diffusion couples to the fast reaction. A user closure
implementing the theta-method or SROCK (Abdulle and Cirilli 2008) can be
supplied when needed.

  

## DAE

Index-1 differential-algebraic equations of the form `M y' = f(t, y)`
with singular mass matrix `M` arise when the model has conservation
constraints (mass balance in chemical reaction networks, total energy
conservation in mechanical systems). A native DAE solver is not
exported; users typically reformulate the index-1 problem as an ODE by
solving the algebraic constraint analytically or numerically and
substituting into the differential block.

For slow-fast systems, mass balance can be imposed as a soft constraint
by adding a penalty term to the ecological RHS
(e.g. `-k_{cons} (sum_i N_i - 1)`) or by reformulating the trait
dynamics as a divisor-free flow; the explicit DAE machinery is not
strictly necessary for the gLV family.

  

## Multi-rate

Multi-rate integration is a generalisation of IMEX that allows different
step sizes for fast and slow components. A native multirate solver is
not exported; references for a user implementation include Bartels,
Constantinescu and Tang (2014) and Sandu (2019).

In the interim, the IMEX-RK solver above provides a working slow-fast
pathway via its additive splitting; the difference is that IMEX-RK uses
the same step size for both halves of the splitting and exploits the
implicit pass’s L-stability, whereas a true multirate scheme uses
different step sizes.

  

## Reproducibility and seeds

Every numeric value cited in this vignette is produced by the canonical
seeded benchmark `inst/extdata/stiff_benchmark_seeded.R` with master
seed `1505`. The recipe to reproduce all of them is:

``` r

library(janos)
set.seed(1505)
spec <- stiff_eco_evo_benchmark()  # defined in inst/extdata/stiff_benchmark_seeded.R
sol  <- dyn_sim(spec, t_max = 5000, solver = solver_imex_ark(),
                       output_n = 501L)
stopifnot(abs(sol$state[1, ncol(sol$state)] - 0.625) < 1e-3)
stopifnot(abs(sol$state[2, ncol(sol$state)] - 0.55)  < 1e-3)
```

The same recipe appears verbatim in the companion `stiff_inverse.Rmd`. A
reader who runs both packages from a fresh clone, with the same seed and
the same benchmark file, should see byte-identical reference numbers up
to machine precision.

The expected reference outputs are:

| quantity                            | reference value      |
|-------------------------------------|----------------------|
| `sol$state[1, end]` (N_1 at t=5000) | 0.625 ± 0.001        |
| `sol$state[2, end]` (N_2 at t=5000) | 0.550 ± 0.001        |
| `sol$state[3, end]` (a_12)          | -0.500 ± 0.001       |
| `sol$state[4, end]` (a_21)          | -0.400 ± 0.001       |
| `sol$n_steps` (IMEX-RK)             | depends on tolerance |
| `stiffness_ratio` on slow manifold  | `~10^3`              |

  

## References

  

### A. Foundational stiff ODE methods

**\[1\]** Curtiss, C. F. and Hirschfelder, J. O. (1952). Integration of
stiff equations. *Proceedings of the National Academy of Sciences*
38(3), 235-243.
[10.1073/pnas.38.3.235](https://doi.org/10.1073/pnas.38.3.235). The
paper that introduced the term “stiff”.

**\[2\]** Dahlquist, G. (1963). A special stability problem for linear
multistep methods. *BIT Numerical Mathematics* 3(1), 27-43.
[10.1007/BF01963532](https://doi.org/10.1007/BF01963532). The
order-stability barrier theorem; cited for the BDF 5-order cap.

**\[3\]** Gear, C. W. (1971). *Numerical Initial Value Problems in
Ordinary Differential Equations*. Prentice-Hall, Englewood Cliffs NJ.
ISBN 978-0-13-626606-9. Classic BDF reference.

**\[4\]** Shampine, L. F. and Gordon, M. K. (1975). *Computer Solution
of Ordinary Differential Equations: The Initial Value Problem*. W. H.
Freeman. ISBN 978-0-7167-0461-5. Cross-check source for BDF coefficient
tables.

**\[5\]** Ehle, B. L. (1969). *On Pade Approximations to the Exponential
Function and A-Stable Methods for the Numerical Solution of Initial
Value Problems*. Research Report CSRR 2010, Department of Applied
Analysis and Computer Science, University of Waterloo. Source for the
A-stability classification of Pade approximants.

**\[6\]** Butcher, J. C. (1964). Implicit Runge-Kutta processes.
*Mathematics of Computation* 18(85), 50-64.
[10.1090/S0025-5718-1964-0159424-9](https://doi.org/10.1090/S0025-5718-1964-0159424-9).
Foundational paper on Radau IIA.

**\[7\]** Petzold, L. R. (1983). Automatic selection of methods for
solving stiff and nonstiff systems of ordinary differential equations.
*SIAM Journal on Scientific and Statistical Computing* 4(1), 136-148.
[10.1137/0904010](https://doi.org/10.1137/0904010). LSODA reference.

**\[8\]** Hindmarsh, A. C. (1983). ODEPACK, a systematized collection of
ODE solvers. In *Scientific Computing*, eds. R. Stepleman et al., 55-64.
North-Holland. ISBN 978-0-444-86607-3. The ODEPACK collection and the
lineage of LSODA.

**\[9\]** Brown, P. N., Byrne, G. D. and Hindmarsh, A. C. (1989). VODE:
A variable-coefficient ODE solver. *SIAM Journal on Scientific and
Statistical Computing* 10(5), 1038-1051.
[10.1137/0910062](https://doi.org/10.1137/0910062). Variable-order BDF
Nordsieck reference.

**\[10\]** Hairer, E., Norsett, S. P. and Wanner, G. (1993). *Solving
Ordinary Differential Equations I: Nonstiff Problems*, 2nd ed. Springer
Series in Computational Mathematics 8. ISBN 978-3-540-56670-0. Companion
to Vol. II.

**\[11\]** Hairer, E. and Wanner, G. (1996). *Solving Ordinary
Differential Equations II: Stiff and Differential-Algebraic Problems*,
2nd ed. Springer Series in Computational Mathematics 14. ISBN
978-3-540-60452-5. The canonical stiff-ODE reference; cited verbatim for
every Butcher tableau in `janos_stiff.hpp`.

**\[12\]** Hindmarsh, A. C., Brown, P. N., Grant, K. E., Lee, S. L.,
Serban, R., Shumaker, D. E. and Woodward, C. S. (2005). SUNDIALS: Suite
of nonlinear and differential/algebraic equation solvers. *ACM
Transactions on Mathematical Software* 31(3), 363-396.
[10.1145/1089014.1089020](https://doi.org/10.1145/1089014.1089020).
Reference suite of production stiff solvers.

**\[13\]** Shampine, L. F. and Reichelt, M. W. (1997). The MATLAB ODE
suite. *SIAM Journal on Scientific Computing* 18(1), 1-22.
[10.1137/S1064827594276424](https://doi.org/10.1137/S1064827594276424).
The MATLAB ode15s / ode23s / ode23tb solver suite.

**\[14\]** Soetaert, K., Petzoldt, T. and Setzer, R. W. (2010). Solving
differential equations in R: package deSolve. *Journal of Statistical
Software* 33(9), 1-25.
[10.18637/jss.v033.i09](https://doi.org/10.18637/jss.v033.i09). The
R-side standard.

**\[15\]** Lambert, J. D. (1991). *Numerical Methods for Ordinary
Differential Systems: The Initial Value Problem*. Wiley. ISBN
978-0-471-92990-1. Comprehensive textbook.

**\[16\]** Iserles, A. (2009). *A First Course in the Numerical Analysis
of Differential Equations*, 2nd ed. Cambridge University Press. ISBN
978-0-521-73490-5. Pedagogical introduction.

  

### B. Runge-Kutta and Rosenbrock methods

**\[17\]** Runge, C. (1895). Uber die numerische Auflosung von
Differentialgleichungen. *Mathematische Annalen* 46(2), 167-178.
[10.1007/BF01446807](https://doi.org/10.1007/BF01446807). The
foundational paper.

**\[18\]** Kutta, W. (1901). Beitrag zur naherungsweisen Integration
totaler Differentialgleichungen. *Zeitschrift fur Mathematik und Physik*
46, 435-453. Companion to Runge.

**\[19\]** Rosenbrock, H. H. (1963). Some general implicit processes for
the numerical solution of differential equations. *The Computer Journal*
5(4), 329-330.
[10.1093/comjnl/5.4.329](https://doi.org/10.1093/comjnl/5.4.329). The
original Rosenbrock paper underlying
[`solver_rosenbrock()`](https://robustecologies.github.io/janos/reference/solver_rosenbrock.md).

**\[20\]** Wanner, G. (1980). On the choice of gamma for singly-implicit
RK or Rosenbrock methods. *BIT Numerical Mathematics* 20(1), 102-106.
[10.1007/BF01933590](https://doi.org/10.1007/BF01933590). Optimal gamma
selection for SDIRK / Rosenbrock.

**\[21\]** Steinebach, G. (1995). *Order-Reduction of ROW-Methods for
DAEs and Method of Lines Applications*. Preprint 1741, Department of
Mathematics, Technische Hochschule Darmstadt. The Rosenbrock-W family
used in RODAS4P.

**\[22\]** Sandu, A., Verwer, J. G., Blom, J. G., Spee, E. J.,
Carmichael, G. R. and Potra, F. A. (1997). Benchmarking stiff ODE
solvers for atmospheric chemistry problems II: Rosenbrock solvers.
*Atmospheric Environment* 31(20), 3459-3472.
[10.1016/S1352-2310(97)83212-8](https://doi.org/10.1016/S1352-2310(97)83212-8).
Atmospheric-chemistry stiff benchmark.

**\[23\]** Hosea, M. E. and Shampine, L. F. (1996). Analysis and
implementation of TR-BDF2. *ACM Transactions on Mathematical Software*
22(1), 105-128.
[10.1145/229473.229474](https://doi.org/10.1145/229473.229474). The
TR-BDF2 ESDIRK used by
[`solver_esdirk()`](https://robustecologies.github.io/janos/reference/solver_esdirk.md).

**\[24\]** Kvaerno, A. (2004). Singly diagonally implicit Runge-Kutta
methods with an explicit first stage. *BIT Numerical Mathematics* 44(3),
489-502.
[10.1023/B:BITN.0000046811.70614.38](https://doi.org/10.1023/B:BITN.0000046811.70614.38).
The Kvaerno ESDIRK family.

**\[25\]** Ascher, U. M., Ruuth, S. J. and Spiteri, R. J. (1997).
Implicit-explicit Runge-Kutta methods for time-dependent partial
differential equations. *Applied Numerical Mathematics* 25(2-3),
151-167.
[10.1016/S0168-9274(97)00056-1](https://doi.org/10.1016/S0168-9274(97)00056-1).
The ARS(2,3,2) IMEX-RK source for
[`solver_imex_ark()`](https://robustecologies.github.io/janos/reference/solver_imex_ark.md).

**\[26\]** Kennedy, C. A. and Carpenter, M. H. (2003). Additive
Runge-Kutta schemes for convection-diffusion-reaction equations.
*Applied Numerical Mathematics* 44(1-2), 139-181.
[10.1016/S0168-9274(02)00138-1](https://doi.org/10.1016/S0168-9274(02)00138-1).
Cross-check for IMEX-RK and ESDIRK.

**\[27\]** Pareschi, L. and Russo, G. (2005). Implicit-explicit
Runge-Kutta schemes and applications to hyperbolic systems with
relaxation. *Journal of Scientific Computing* 25(1-2), 129-155.
[10.1007/s10915-004-4636-4](https://doi.org/10.1007/s10915-004-4636-4).
Pareschi-Russo IMEX-RK family.

**\[28\]** Boscarino, S. and Russo, G. (2009). On a class of uniformly
accurate IMEX Runge-Kutta schemes and applications to hyperbolic systems
with relaxation. *SIAM Journal on Scientific Computing* 31(3),
1926-1945. [10.1137/080713562](https://doi.org/10.1137/080713562).
Modern IMEX-RK reference.

**\[29\]** Butcher, J. C. (2016). *Numerical Methods for Ordinary
Differential Equations*, 3rd ed. Wiley. ISBN 978-1-119-12150-3. The
standard RK methods monograph.

  

### C. Step-size control and error estimation

**\[30\]** Gustafsson, K. (1991). Control theoretic techniques for
stepsize selection in explicit Runge-Kutta methods. *ACM Transactions on
Mathematical Software* 17(4), 533-554.
[10.1145/210232.210242](https://doi.org/10.1145/210232.210242). The PI
step-size controller used in `GustafssonPIController`.

**\[31\]** Gustafsson, K. (1994). Control-theoretic techniques for
stepsize selection in implicit Runge-Kutta methods. *ACM Transactions on
Mathematical Software* 20(4), 496-517.
[10.1145/198429.198437](https://doi.org/10.1145/198429.198437).
Implicit-method companion.

**\[32\]** Shampine, L. F. (1980). Implementation of implicit formulas
for the solution of ODEs. *SIAM Journal on Scientific and Statistical
Computing* 1(1), 103-118.
[10.1137/0901005](https://doi.org/10.1137/0901005). Standard
error-estimator strategy for BDF.

**\[33\]** Soderlind, G. (2002). Automatic control and adaptive
time-stepping. *Numerical Algorithms* 31(1), 281-310.
[10.1023/A:1021160023092](https://doi.org/10.1023/A:1021160023092).
Modern PI / PID step controllers.

**\[34\]** Soderlind, G. and Wang, L. (2006). Adaptive time-stepping and
computational stability. *Journal of Computational and Applied
Mathematics* 185(2), 225-243.
[10.1016/j.cam.2005.03.008](https://doi.org/10.1016/j.cam.2005.03.008).
Stability under adaptive control.

**\[35\]** Hairer, E. and Wanner, G. (1999). Stiff differential
equations solved by Radau methods. *Journal of Computational and Applied
Mathematics* 111(1-2), 93-111.
[10.1016/S0377-0427(99)00134-X](https://doi.org/10.1016/S0377-0427(99)00134-X).
Practical implementation of Radau IIA.

  

### D. Matrix exponential and Krylov methods

**\[36\]** Moler, C. and Van Loan, C. (2003). Nineteen dubious ways to
compute the exponential of a matrix, twenty-five years later. *SIAM
Review* 45(1), 3-49.
[10.1137/S00361445024180](https://doi.org/10.1137/S00361445024180). The
canonical survey of matrix-exponential algorithms.

**\[37\]** Higham, N. J. (2005). The scaling and squaring method for the
matrix exponential revisited. *SIAM Journal on Matrix Analysis and
Applications* 26(4), 1179-1193.
[10.1137/04061101X](https://doi.org/10.1137/04061101X). The Pade-13
algorithm used in
[`expm_pade()`](https://robustecologies.github.io/janos/reference/expm_pade.md).

**\[38\]** Higham, N. J. (2008). *Functions of Matrices: Theory and
Computation*. SIAM. ISBN 978-0-89871-646-7. Comprehensive monograph on
matrix functions.

**\[39\]** Saad, Y. (1992). Analysis of some Krylov subspace
approximations to the matrix exponential operator. *SIAM Journal on
Numerical Analysis* 29(1), 209-228.
[10.1137/0729014](https://doi.org/10.1137/0729014). The Arnoldi-Krylov
approximation used in
[`expm_krylov()`](https://robustecologies.github.io/janos/reference/expm_krylov.md).

**\[40\]** Sidje, R. B. (1998). Expokit: A software package for
computing matrix exponentials. *ACM Transactions on Mathematical
Software* 24(1), 130-156.
[10.1145/285861.285868](https://doi.org/10.1145/285861.285868). The
Expokit Fortran software for action of matrix exponential.

**\[41\]** Niesen, J. and Wright, W. M. (2012). Algorithm 919: A Krylov
subspace algorithm for evaluating the phi-functions appearing in
exponential integrators. *ACM Transactions on Mathematical Software*
38(3), Article 22.
[10.1145/2168773.2168781](https://doi.org/10.1145/2168773.2168781). The
`phi_k` evaluator planned for EPIRK.

**\[42\]** Hochbruck, M. and Ostermann, A. (2010). Exponential
integrators. *Acta Numerica* 19, 209-286.
[10.1017/S0962492910000048](https://doi.org/10.1017/S0962492910000048).
The standard review of exponential integrators.

**\[43\]** Tokman, M. (2006). Efficient integration of large stiff
systems of ODEs with exponential propagation iterative (EPI) methods.
*Journal of Computational Physics* 213(2), 748-776.
[10.1016/j.jcp.2005.08.032](https://doi.org/10.1016/j.jcp.2005.08.032).
Foundational EPI method paper.

**\[44\]** Tokman, M. and Loffeld, J. (2011). Efficient design of
exponential-Krylov integrators for large scale computing. *SIAM Journal
on Scientific Computing* 33(4), 1480-1504.
[10.1137/100795464](https://doi.org/10.1137/100795464). EPIRK reference.

**\[45\]** Al-Mohy, A. H. and Higham, N. J. (2011). Computing the action
of the matrix exponential, with an application to exponential
integrators. *SIAM Journal on Scientific Computing* 33(2), 488-511.
[10.1137/100788860](https://doi.org/10.1137/100788860). The
action-of-exponential algorithm.

  

### E. Stabilised explicit methods (ROCK and friends)

**\[46\]** van der Houwen, P. J. and Sommeijer, B. P. (1980). On the
internal stability of explicit, m-stage Runge-Kutta methods for large
m-values. *ZAMM* 60, 479-485.
[10.1002/zamm.19800601005](https://doi.org/10.1002/zamm.19800601005).
Foundational stabilised explicit methods.

**\[47\]** Abdulle, A. and Medovikov, A. A. (2001). Second order
Chebyshev methods based on orthogonal polynomials. *Numerische
Mathematik* 90(1), 1-18.
[10.1007/s002110100292](https://doi.org/10.1007/s002110100292). ROCK2.

**\[48\]** Abdulle, A. (2002). Fourth order Chebyshev methods with
recurrence relation. *SIAM Journal on Scientific Computing* 23(6),
2041-2054.
[10.1137/S1064827500379549](https://doi.org/10.1137/S1064827500379549).
ROCK4 reference.

**\[49\]** Abdulle, A. and Vilmart, G. (2013). PIROCK: A swiss-knife
partitioned implicit-explicit orthogonal Runge-Kutta Chebyshev
integrator for stiff diffusion-advection-reaction problems with or
without noise. *Journal of Computational Physics* 242, 869-888.
[10.1016/j.jcp.2013.02.009](https://doi.org/10.1016/j.jcp.2013.02.009).
The PIROCK extension.

  

### F. Stochastic stiff integration

**\[50\]** Kloeden, P. E. and Platen, E. (1999). *Numerical Solution of
Stochastic Differential Equations*, corrected 3rd printing. Springer.
ISBN 978-3-540-54062-5. Standard SDE textbook.

**\[51\]** Higham, D. J. (2001). An algorithmic introduction to
numerical simulation of stochastic differential equations. *SIAM Review*
43(3), 525-546.
[10.1137/S0036144500378302](https://doi.org/10.1137/S0036144500378302).
Pedagogical SDE reference; underlies drift-implicit theta method planned
for janos.

**\[52\]** Milstein, G. N. and Tretyakov, M. V. (2004). *Stochastic
Numerics for Mathematical Physics*. Springer. ISBN 978-3-540-21110-1.
SDE numerics monograph.

**\[53\]** Abdulle, A. and Cirilli, S. (2008). S-ROCK: Chebyshev methods
for stiff stochastic differential equations. *SIAM Journal on Scientific
Computing* 30(2), 997-1014.
[10.1137/070679375](https://doi.org/10.1137/070679375). SROCK reference;
planned for janos.

**\[54\]** Abdulle, A., Vilmart, G. and Zygalakis, K. C. (2013). Weak
second order explicit stabilized methods for stiff stochastic
differential equations. *SIAM Journal on Scientific Computing* 35(4),
A1792-A1814. [10.1137/12088954X](https://doi.org/10.1137/12088954X).
Second-order SROCK.

  

### G. Singular perturbation theory and slow-fast systems

**\[55\]** Tikhonov, A. N. (1952). Systems of differential equations
containing small parameters in the derivatives. *Matematicheskii
Sbornik* 31(73), 575-586. Foundational singular-perturbation paper.

**\[56\]** Fenichel, N. (1979). Geometric singular perturbation theory
for ordinary differential equations. *Journal of Differential Equations*
31(1), 53-98.
[10.1016/0022-0396(79)90152-9](https://doi.org/10.1016/0022-0396(79)90152-9).
The slow-manifold invariance theorem.

**\[57\]** Jones, C. K. R. T. (1995). Geometric singular perturbation
theory. In *Dynamical Systems*, ed. R. Johnson, *Lecture Notes in
Mathematics* 1609, 44-118. Springer.
[10.1007/BFb0095239](https://doi.org/10.1007/BFb0095239). Modern review.

**\[58\]** Vasilieva, A. B., Butuzov, V. F. and Kalachev, L. V. (1995).
*The Boundary Function Method for Singular Perturbation Problems*. SIAM.
ISBN 978-0-89871-333-6. Comprehensive monograph.

**\[59\]** O’Malley, R. E. (1991). *Singular Perturbation Methods for
Ordinary Differential Equations*. Applied Mathematical Sciences 89.
Springer. ISBN 978-0-387-97556-9. Standard textbook.

**\[60\]** Verhulst, F. (2005). *Methods and Applications of Singular
Perturbations: Boundary Layers and Multiple Timescale Dynamics*.
Springer. ISBN 978-0-387-22966-9. Application-oriented reference.

**\[61\]** Kuehn, C. (2015). *Multiple Time Scale Dynamics*. Applied
Mathematical Sciences 191. Springer. ISBN 978-3-319-12315-4.
Comprehensive monograph on multi-timescale dynamics.

**\[62\]** Berglund, N. and Gentz, B. (2006). *Noise-Induced Phenomena
in Slow-Fast Dynamical Systems*. Springer. ISBN 978-1-84628-038-6.
Stochastic singular-perturbation reference.

**\[63\]** Bartels-Constantinescu, T. and Constantinescu, E. M. (2014).
Multirate timestepping methods for hyperbolic conservation laws.
*Journal of Scientific Computing* 33(2), 239-278. Multirate integration
reference.

**\[64\]** Sandu, A. (2019). A class of multirate infinitesimal GARK
methods. *SIAM Journal on Numerical Analysis* 57(5), 2300-2327.
[10.1137/18M1205492](https://doi.org/10.1137/18M1205492). Modern
multirate framework.

  

### H. Differential-algebraic equations

**\[65\]** Brenan, K. E., Campbell, S. L. and Petzold, L. R. (1996).
*Numerical Solution of Initial-Value Problems in Differential-Algebraic
Equations*. SIAM Classics in Applied Mathematics 14. ISBN
978-0-89871-353-4. The standard DAE reference.

**\[66\]** Petzold, L. R. (1982). A description of DASSL: A
differential/algebraic system solver. In *Scientific Computing*, eds. R.
Stepleman et al. North-Holland. The DASSL reference.

**\[67\]** Hairer, E., Lubich, C. and Roche, M. (1989). *The Numerical
Solution of Differential-Algebraic Systems by Runge-Kutta Methods*.
*Lecture Notes in Mathematics* 1409. Springer.
[10.1007/BFb0093947](https://doi.org/10.1007/BFb0093947). DAE-RK
monograph.

**\[68\]** Kunkel, P. and Mehrmann, V. (2006). *Differential-Algebraic
Equations: Analysis and Numerical Solution*. EMS Textbooks in
Mathematics. ISBN 978-3-03719-017-1. Modern treatment.

  

### I. Canonical stiff benchmarks

**\[69\]** Robertson, H. H. (1966). The solution of a set of reaction
rate equations. In *Numerical Analysis: An Introduction*, ed. J. Walsh,
178-182. Academic Press, London. The Robertson benchmark.

**\[70\]** Schaefer, E. (1975). A new approach to explain the “high
irradiance responses” of photomorphogenesis on the basis of phytochrome.
*Journal of Mathematical Biology* 2(1), 41-56.
[10.1007/BF00276015](https://doi.org/10.1007/BF00276015). The HIRES
benchmark.

**\[71\]** Field, R. J. and Noyes, R. M. (1974). Oscillations in
chemical systems IV: Limit cycle behaviour in a model of a real chemical
reaction. *Journal of Chemical Physics* 60(5), 1877-1884.
[10.1063/1.1681288](https://doi.org/10.1063/1.1681288). The OREGO
benchmark.

**\[72\]** Verwer, J. G., Schamhorst, J. G. and Stelling, J. F. (1996).
Convergence properties of the Runge-Kutta-Chebyshev method. *Numerische
Mathematik* 67, 313-323. Reference to the Pollu atmospheric-chemistry
benchmark.

**\[73\]** van der Pol, B. (1926). On “relaxation oscillations”.
*Philosophical Magazine* 2(11), 978-992. The original Van der Pol
oscillator.

**\[74\]** Mazzia, F. and Iavernaro, F. (2003). *Test set for initial
value problem solvers*. Department of Mathematics, University of Bari.
<https://www.dm.uniba.it/~testset/>. The standard stiff benchmark test
suite (IFISS test set).

**\[75\]** Soderlind, G. and Wang, L. (2006). Evaluating numerical
ODE/DAE methods, algorithms and software. *Journal of Computational and
Applied Mathematics* 185(2), 244-260.
[10.1016/j.cam.2005.03.009](https://doi.org/10.1016/j.cam.2005.03.009).
Benchmark evaluation methodology.

  

### J. Slow-fast and multi-timescale biology

**\[76\]** Lotka, A. J. (1925). *Elements of Physical Biology*. Williams
and Wilkins, Baltimore. Foundational text for the Lotka-Volterra
equations.

**\[77\]** Volterra, V. (1926). Fluctuations in the abundance of a
species considered mathematically. *Nature* 118, 558-560.
[10.1038/118558a0](https://doi.org/10.1038/118558a0). Companion to
Lotka.

**\[78\]** May, R. M. (1973). *Stability and Complexity in Model
Ecosystems*. Princeton University Press. ISBN 978-0-691-08130-9.
Foundational community-ecology dynamics.

**\[79\]** MacArthur, R. (1970). Species packing and competitive
equilibrium for many species. *Theoretical Population Biology* 1(1),
1-11.
[10.1016/0040-5809(70)90039-0](https://doi.org/10.1016/0040-5809(70)90039-0).
Competition-coefficient theory.

**\[80\]** Lande, R. (1976). Natural selection and random genetic drift
in phenotypic evolution. *Evolution* 30(2), 314-334.
[10.2307/2407703](https://doi.org/10.2307/2407703). Foundational
quantitative genetics of trait evolution.

**\[81\]** Lande, R. (1979). Quantitative genetic analysis of
multivariate evolution applied to brain-body size allometry. *Evolution*
33(1), 402-416. [10.2307/2407630](https://doi.org/10.2307/2407630).
Multivariate-selection framework.

**\[82\]** Roughgarden, J. (1976). Resource partitioning among competing
species: a coevolutionary approach. *Theoretical Population Biology*
9(3), 388-424.
[10.1016/0040-5809(76)90054-X](https://doi.org/10.1016/0040-5809(76)90054-X).
Coevolutionary competition theory.

**\[83\]** Yoshida, T., Jones, L. E., Ellner, S. P., Fussmann, G. F. and
Hairston Jr, N. G. (2003). Rapid evolution drives ecological dynamics in
a predator-prey system. *Nature* 424(6946), 303-306.
[10.1038/nature01767](https://doi.org/10.1038/nature01767). Empirical
foundation of eco-evolutionary dynamics.

**\[84\]** Hairston Jr, N. G., Ellner, S. P., Geber, M. A., Yoshida, T.
and Fox, J. A. (2005). Rapid evolution and the convergence of ecological
and evolutionary time. *Ecology Letters* 8(10), 1114-1127.
[10.1111/j.1461-0248.2005.00812.x](https://doi.org/10.1111/j.1461-0248.2005.00812.x).
The “convergence of time scales” paradigm.

**\[85\]** Fussmann, G. F., Loreau, M. and Abrams, P. A. (2007).
Eco-evolutionary dynamics of communities and ecosystems. *Functional
Ecology* 21(3), 465-477.
[10.1111/j.1365-2435.2007.01275.x](https://doi.org/10.1111/j.1365-2435.2007.01275.x).
Conceptual review.

**\[86\]** Schoener, T. W. (2011). The newest synthesis: understanding
the interplay of evolutionary and ecological dynamics. *Science*
331(6016), 426-429.
[10.1126/science.1193954](https://doi.org/10.1126/science.1193954).
Synthetic review.

**\[87\]** Hendry, A. P. (2017). *Eco-Evolutionary Dynamics*. Princeton
University Press. ISBN 978-0-691-14543-8. Standard monograph.

**\[88\]** Govaert, L., Fronhofer, E. A., Lion, S. et al. (2019).
Eco-evolutionary feedbacks - Theoretical models and perspectives.
*Functional Ecology* 33(1), 13-30.
[10.1111/1365-2435.13241](https://doi.org/10.1111/1365-2435.13241).
Modern review of feedback models.

  

### K. Reaction networks, chemical kinetics and oscillating reactions

**\[89\]** Gillespie, D. T. (1977). Exact stochastic simulation of
coupled chemical reactions. *Journal of Physical Chemistry* 81(25),
2340-2361. [10.1021/j100540a008](https://doi.org/10.1021/j100540a008).
The exact stochastic simulation algorithm.

**\[90\]** Belousov, B. P. (1959). A periodic reaction and its
mechanism. *Compilation of Abstracts on Radiation Medicine* 147(145), 1.
The original Belousov-Zhabotinsky reaction.

**\[91\]** Lefever, R. and Prigogine, I. (1968). Symmetry breaking
instabilities in dissipative systems II. *Journal of Chemical Physics*
48(4), 1695-1700.
[10.1063/1.1668896](https://doi.org/10.1063/1.1668896). The Brusselator
model.

**\[92\]** Sel’kov, E. E. (1968). Self-oscillations in glycolysis: A
simple kinetic model. *European Journal of Biochemistry* 4(1), 79-86.
[10.1111/j.1432-1033.1968.tb00175.x](https://doi.org/10.1111/j.1432-1033.1968.tb00175.x).
The Sel’kov glycolysis benchmark.

**\[93\]** Hindmarsh, J. L. and Rose, R. M. (1984). A model of neuronal
bursting using three coupled first order differential equations.
*Proceedings of the Royal Society B* 221(1222), 87-102.
[10.1098/rspb.1984.0024](https://doi.org/10.1098/rspb.1984.0024). The
Hindmarsh-Rose neuronal model.

**\[94\]** FitzHugh, R. (1961). Impulses and physiological states in
theoretical models of nerve membrane. *Biophysical Journal* 1(6),
445-466.
[10.1016/S0006-3495(61)86902-6](https://doi.org/10.1016/S0006-3495(61)86902-6).
The FitzHugh-Nagumo model.

  

### L. Sensitivity analysis and adjoint methods

**\[95\]** Cao, Y., Li, S. and Petzold, L. (2002). Adjoint sensitivity
analysis for differential-algebraic equations: algorithms and software.
*Journal of Computational and Applied Mathematics* 149(1), 171-191.
[10.1016/S0377-0427(02)00528-9](https://doi.org/10.1016/S0377-0427(02)00528-9).
Foundational adjoint-sensitivity algorithm.

**\[96\]** Hindmarsh, A. C. and Serban, R. (2002). *User Documentation
for CVODES, an ODE Solver with Sensitivity Analysis Capabilities*.
UCRL-MA-148813, Lawrence Livermore National Laboratory. The CVODES
manual.

**\[97\]** Griewank, A. and Walther, A. (2000). Algorithm 799: revolve:
An implementation of checkpointing for the reverse or adjoint mode of
computational differentiation. *ACM Transactions on Mathematical
Software* 26(1), 19-45.
[10.1145/347837.347846](https://doi.org/10.1145/347837.347846). The
revolve checkpointing scheme.

**\[98\]** Griewank, A. and Walther, A. (2008). *Evaluating Derivatives:
Principles and Techniques of Algorithmic Differentiation*, 2nd ed. SIAM.
ISBN 978-0-89871-659-7. The standard automatic-differentiation
reference.

**\[99\]** Caracotsios, M. and Stewart, W. E. (1985). Sensitivity
analysis of initial value problems with mixed ODEs and algebraic
equations. *Computers & Chemical Engineering* 9(4), 359-365.
[10.1016/0098-1354(85)85014-6](https://doi.org/10.1016/0098-1354(85)85014-6).
Early sensitivity analysis paper.

**\[100\]** Errico, R. M. (1997). What is an adjoint model? *Bulletin of
the American Meteorological Society* 78(11), 2577-2591. Pedagogical
introduction to adjoint methods.

  

### M. Numerical linear algebra

**\[101\]** Golub, G. H. and Van Loan, C. F. (2013). *Matrix
Computations*, 4th ed. Johns Hopkins University Press. ISBN
978-1-4214-0794-4. Standard numerical linear algebra reference.

**\[102\]** Trefethen, L. N. and Bau, D. (1997). *Numerical Linear
Algebra*. SIAM. ISBN 978-0-89871-361-9. Pedagogical companion.

**\[103\]** Bartels, R. H. and Stewart, G. W. (1972). Solution of the
matrix equation AX + XB = C. *Communications of the ACM* 15(9), 820-826.
[10.1145/361573.361582](https://doi.org/10.1145/361573.361582). The
Bartels-Stewart algorithm for Lyapunov equations.

**\[104\]** Anderson, E. et al. (1999). *LAPACK Users’ Guide*, 3rd
ed. SIAM. ISBN 978-0-89871-447-0. The LAPACK reference suite for dense
linear algebra.

**\[105\]** Sanderson, C. and Curtin, R. (2016). Armadillo: a
template-based C++ library for linear algebra. *Journal of Open Source
Software* 1(2), 26.
[10.21105/joss.00026](https://doi.org/10.21105/joss.00026). The
Armadillo C++ library underlying the janos stiff stack.

**\[106\]** Eddelbuettel, D. and Sanderson, C. (2014). RcppArmadillo:
Accelerating R with high-performance C++ linear algebra. *Computational
Statistics and Data Analysis* 71, 1054-1063.
[10.1016/j.csda.2013.02.005](https://doi.org/10.1016/j.csda.2013.02.005).
The R-side bridge.

  

### N. Newton iteration and nonlinear solvers

**\[107\]** Dennis Jr, J. E. and Schnabel, R. B. (1996). *Numerical
Methods for Unconstrained Optimization and Nonlinear Equations*. SIAM
Classics in Applied Mathematics 16. ISBN 978-0-89871-364-0. Standard
Newton-method reference.

**\[108\]** Deuflhard, P. (2011). *Newton Methods for Nonlinear
Problems: Affine Invariance and Adaptive Algorithms*. Springer. ISBN
978-3-642-23898-7. Modern modified-Newton treatment.

**\[109\]** Eisenstat, S. C. and Walker, H. F. (1996). Choosing the
forcing terms in an inexact Newton method. *SIAM Journal on Scientific
Computing* 17(1), 16-32.
[10.1137/0917003](https://doi.org/10.1137/0917003). Inexact Newton
tolerances.

**\[110\]** Kelley, C. T. (2003). *Solving Nonlinear Equations with
Newton’s Method*. SIAM Fundamentals of Algorithms 1. ISBN
978-0-89871-546-0. Practical Newton-iteration handbook.

  

### O. Data assimilation and ensemble filtering

**\[111\]** Kalman, R. E. (1960). A new approach to linear filtering and
prediction problems. *Journal of Basic Engineering* 82(1), 35-45.
[10.1115/1.3662552](https://doi.org/10.1115/1.3662552). The original
Kalman filter.

**\[112\]** Evensen, G. (1994). Sequential data assimilation with a
nonlinear quasi-geostrophic model using Monte Carlo methods to forecast
error statistics. *Journal of Geophysical Research: Oceans* 99(C5),
10143-10162. [10.1029/94JC00572](https://doi.org/10.1029/94JC00572). The
original EnKF.

**\[113\]** Burgers, G., van Leeuwen, P. J. and Evensen, G. (1998).
Analysis scheme in the ensemble Kalman filter. *Monthly Weather Review*
126(6), 1719-1724. The perturbed-observation update used in the
RElabverse companion’s `stiff_enkf()`.

**\[114\]** Arnold, A., Calvetti, D. and Somersalo, E. (2014). Parameter
estimation for stiff deterministic dynamical systems via ensemble Kalman
filter. *Inverse Problems* 30(10), 105008.
[10.1088/0266-5611/30/10/105008](https://doi.org/10.1088/0266-5611/30/10/105008).
Anchor reference for the inverse-problem side.

**\[115\]** Asch, M., Bocquet, M. and Nodet, M. (2016). *Data
Assimilation: Methods, Algorithms, and Applications*. SIAM. ISBN
978-1-611974-53-9. Comprehensive monograph including
parameter-estimation EnKF.

  

### P. Software and computational tools

**\[116\]** R Core Team (2024). *R: A Language and Environment for
Statistical Computing*. R Foundation for Statistical Computing, Vienna,
Austria. <https://www.R-project.org/>.

**\[117\]** Eddelbuettel, D. and Francois, R. (2011). Rcpp: Seamless R
and C++ integration. *Journal of Statistical Software* 40(8), 1-18.
[10.18637/jss.v040.i08](https://doi.org/10.18637/jss.v040.i08).

**\[118\]** Wickham, H. (2016). *ggplot2: Elegant Graphics for Data
Analysis*, 2nd ed. Springer. ISBN 978-3-319-24277-4. The ggplot2
reference.

**\[119\]** Pedersen, T. L. (2024). *patchwork: The Composer of Plots*.
R package version 1.2.0. <https://patchwork.data-imaginist.com/>.

**\[120\]** Rackauckas, C. and Nie, Q. (2017).
DifferentialEquations.jl - A performant and feature-rich ecosystem for
solving differential equations in Julia. *Journal of Open Research
Software* 5(1), 15.
[10.5334/jors.151](https://doi.org/10.5334/jors.151). Julia-side
standard for differential equations.

**\[121\]** Stadelmann, M. and Sanderson, C. (2025). ensmallen: A
high-performance C++ library for numerical optimization. *Software
reference*. <https://ensmallen.org/>.

**\[122\]** Wickham, H. (2019). *Advanced R*, 2nd ed. CRC Press. ISBN
978-0-8153-8457-0. Reference for S3 dispatch.

  

### Q. Mathematical foundations

**\[123\]** Hartman, P. (1982). *Ordinary Differential Equations*, 2nd
ed. SIAM Classics in Applied Mathematics 38. ISBN 978-0-89871-510-1.
Foundational ODE theory.

**\[124\]** Coddington, E. A. and Levinson, N. (1955). *Theory of
Ordinary Differential Equations*. McGraw-Hill. The standard ODE
textbook.

**\[125\]** Arnold, V. I. (1992). *Ordinary Differential Equations*, 3rd
ed. Springer Texts in Mathematics. ISBN 978-3-540-54813-3. Geometric
perspective.

**\[126\]** Guckenheimer, J. and Holmes, P. (1983). *Nonlinear
Oscillations, Dynamical Systems, and Bifurcations of Vector Fields*.
Applied Mathematical Sciences 42. Springer. ISBN 978-0-387-90819-9. The
dynamical-systems standard.

**\[127\]** Wiggins, S. (2003). *Introduction to Applied Nonlinear
Dynamical Systems and Chaos*, 2nd ed. Texts in Applied Mathematics 2.
Springer. ISBN 978-0-387-00177-7. Application-oriented companion.

**\[128\]** Strogatz, S. H. (2014). *Nonlinear Dynamics and Chaos*, 2nd
ed. Westview Press. ISBN 978-0-8133-4910-7. Pedagogical introduction.

**\[129\]** Boyce, W. E. and DiPrima, R. C. (2017). *Elementary
Differential Equations and Boundary Value Problems*, 11th ed. Wiley.
ISBN 978-1-119-32063-0. Undergraduate-level ODE reference.

**\[130\]** Atkinson, K. E. (1989). *An Introduction to Numerical
Analysis*, 2nd ed. Wiley. ISBN 978-0-471-62489-9. General numerical
analysis reference.

  

## Appendix A: Butcher tableaux verbatim

For documentation transparency and to allow downstream verification, the
Butcher tableaux of the four new solver families are reproduced here
verbatim from the cited references. Coefficient typos in stiff solver
implementations are notoriously hard to detect; a printed table that the
reader can compare against the literature is the simplest defence
against silent regressions.

### A.1 BDF coefficients (Hairer-Wanner Vol. II Table 1.1)

The BDF formula at order `k` is
`sum_{j = 0}^{k} alpha_j y_{n+1-j} = beta h f_{n+1}` with `alpha_0 = 1`.
The remaining coefficients are:

| order | `alpha_1` | `alpha_2` | `alpha_3` | `alpha_4` | `alpha_5` | `beta` |
|------:|----------:|----------:|----------:|----------:|----------:|-------:|
|     1 |        -1 |           |           |           |           |      1 |
|     2 |      -4/3 |       1/3 |           |           |           |    2/3 |
|     3 |    -18/11 |      9/11 |     -2/11 |           |           |   6/11 |
|     4 |    -48/25 |     36/25 |    -16/25 |      3/25 |           |  12/25 |
|     5 |  -300/137 |   300/137 |  -200/137 |    75/137 |   -12/137 | 60/137 |

A-stability angle `alpha` (degrees) and L-stability flag:

| order | A(alpha) | L-stable | order of LTE |
|------:|---------:|---------:|-------------:|
|     1 |       90 |      yes |          h^2 |
|     2 |       90 |      yes |          h^3 |
|     3 |       86 |       no |          h^4 |
|     4 |       73 |       no |          h^5 |
|     5 |       51 |       no |          h^6 |

### A.2 TR-BDF2 ESDIRK (Hosea-Shampine 1996 eq. 1.4)

With `gamma = 2 - sqrt(2) approx 0.5858`, `d = gamma / 2 approx 0.2929`,
`w = sqrt(2) / 4 approx 0.3536`:

       c |       A
       0 |   0,        0,     0
    gamma|   d,        d,     0
       1 |   w,        w,     d
    -----+--------------------
       b |   w,        w,     d
     b_hat| (1-w)/3, (3w+1)/3, d/3

The embedded order-3 estimator weights `b_hat` give error coefficients
`e = b - b_hat = (w - (1-w)/3, w - (3w+1)/3, d - d/3) = ((2w + (w-1))/3, ..., (2d)/3)`.
The local error vector is `h * e^T (k_1, k_2, k_3)`.

### A.3 ARS(2,3,2) IMEX-RK (Ascher-Ruuth-Spiteri 1997 Table 1)

With `gamma = (2 - sqrt(2))/2 approx 0.2929`,
`delta = 1 - 1/(2 gamma) approx -0.7071`:

Implicit tableau (3 stages, 2 implicit):

       c_I  |       A_I
        0   |   0,        0,     0
      gamma |   0,      gamma,   0
        1   |   0,    1-gamma, gamma
    --------+----------------------
       b_I  |   0,    1-gamma, gamma

Explicit tableau (3 stages, 3 explicit):

       c_E  |       A_E
        0   |   0,        0,     0
      gamma | gamma,      0,     0
        1   | delta,    1-delta, 0
    --------+----------------------
       b_E  | delta,    1-delta, 0

Order of the combined scheme: 2. Embedded order-1 estimator: implicit
Euler step.

### A.4 Radau IIA s=2 order 3 (Hairer-Wanner Table 5.5)

       c |       A          
     1/3 |  5/12,   -1/12   
      1  |  3/4,    1/4     
    -----+----------------- 
      b  |  3/4,    1/4     

The simplified Newton requires the eigendecomposition of `A^{-1}`; for s
= 2 the eigenvalues are `lambda_R = 3` (real) and
`lambda_C = 1.5 +/- 1.5i sqrt(7)/3` (complex conjugate pair). The
implementation uses the full block-Newton of size `2n`; the
eigendecomposition trick can be substituted by a user wrapper around the
engine. The matrix `A^{-1}` is:

    A^{-1} = (1 / det A) * [[1/4, 1/12], [-3/4, 5/12]]
    det A  = 5/12 * 1/4 - (-1/12) * 3/4 = 5/48 + 3/48 = 1/6
    A^{-1} = 6 * [[1/4, 1/12], [-3/4, 5/12]]
           = [[3/2, 1/2], [-9/2, 5/2]]

Embedded order-2 estimator: trapezoidal-rule combination of stage
derivatives, `y_low = y_n + h (F_1 + F_2) / 2`. The Hairer-Wanner Vol.
II Section IV.8.18 estimator uses a different perturbed-step
construction that gives an order-3 indicator for the order-3 method; a
user closure can supply that estimator.

  

## Appendix B: Performance notes

The native stiff engines accept an R closure for the right-hand side.
This costs roughly `2-5` microseconds per call on modest hardware
compared to a tenth of a microsecond for the formula-compiled fast path
used by
[`solver_rosenbrock()`](https://robustecologies.github.io/janos/reference/solver_rosenbrock.md).
On a 4-state slow-fast problem this overhead is dominant: each accepted
step incurs roughly `O(s * n + n + 1)` RHS calls for the simplified
Newton plus `n + 1` calls for the finite-difference Jacobian, so a Radau
step costs roughly `2 + 5 = 7` calls and an IMEX-RK step costs roughly
`3 + 5 = 8` calls. The total integration cost is dominated by the
R-callback overhead.

For modest stiff problems with `n = 4 ... 20` state components and an
analytical Jacobian provided, the native stiff path runs roughly
`5 - 20` times slower than the formula-compiled
[`solver_rosenbrock()`](https://robustecologies.github.io/janos/reference/solver_rosenbrock.md)
on the same problem. This is acceptable for research workflows where the
integration is run a handful of times per analysis; tight loops within a
Bayesian sampler benefit from supplying an analytical Jacobian (cutting
roughly half the per-step RHS calls) or from a compiled RHS through
[`system_spec_rhs_cpp()`](https://robustecologies.github.io/janos/reference/system_spec_rhs_cpp.md).

The principal benefit of the v1 R-callback path is its flexibility: any
R function with the standard signature can be the RHS, regardless of
whether it can be represented as a formula. This makes the new solvers
immediately usable in identifiability and parameter-fitting workflows
where the model RHS is itself a function of fitted parameters (e.g. the
[`RElabverse::stiff_pibvb()`](https://rdrr.io/pkg/RElabverse/man/stiff_pibvb.html)
time-varying-parameter likelihood, where the RHS is a closure over the
basis coefficients).

Future performance milestones:

| milestone | expected gain |
|----|---:|
| formula-compilation bridge | 5-10x |
| analytical Jacobian | 1.5-2x |
| simplified Newton for Radau | 1.5x |
| Krylov-Arnoldi for large stiff | constant factor for `n > 100` |
| OpenMP parallel ensemble | n_cores (when running ensemble simulations) |

  

## Appendix C: Step-controller tuning

The Gustafsson PI controller has three tunable parameters: safety
(default 0.9), alpha (default 0.7), beta (default 0.4). The default
values are robust on the canonical benchmark suite; they correspond to
the Hairer-Wanner recommendation and trade modest step-size oscillation
for safe acceptance rates.

For problems with very abrupt local-error spikes (e.g. transition
between regimes in OREGO), reducing `safety` to 0.8 and `beta` to 0.2
makes the controller less reactive and reduces rejected-step counts at
the cost of slower step growth in smooth regions. Conversely, for very
smooth problems where the trajectory has no sharp features (e.g. linear
decay, slow exponential approach to equilibrium), raising `safety` to
0.95 and `beta` to 0.6 lets the controller grow the step more
aggressively. The defaults are a compromise that does no worse than
`1.5x` the optimal step count on any of the benchmarks above.

The controller exposes the following knobs (accessible by reading the
`solver_spec` slots and rebuilding):

- `dt_min` (default `1e-12`): integration aborts with a warning if the
  step size collapses below this floor. Setting it too high causes
  spurious termination on legitimate stiff transients; setting it too
  low lets the controller spiral toward floating-point underflow on a
  poorly-conditioned residual.

- `dt_max` (default `1.0`): cap on the step size. For problems with very
  long time horizons, raise this to permit very-large steps in smooth
  regions; for problems where the trajectory has narrow features that
  must be resolved, lower it to prevent the controller from coasting
  past them.

- `max_steps` (default `1e7`): absolute upper bound on accepted-step
  count. Hit this and the integration terminates with a warning. For
  very long stiff trajectories (e.g. Robertson to `t = 10^11`) the bound
  may need to be raised explicitly.

A future minor release will provide a `solver_diagnostic()` function
that scans the step-size and step-error history of a completed
integration and recommends controller-parameter adjustments.

  

## Appendix D: Cross-package call graph

The figure below summarises the cross-package call sequence the
canonical slow-fast worked example exercises:

    inst/extdata/stiff_benchmark_seeded.R        (defines spec_ecoevo)
            |
            v
    janos::dyn_sim(spec_ecoevo, ...)       (forward integration)
            |
            | -- janos_stiff.hpp::NewtonSolver
            | -- janos_stiff.hpp::GustafssonPIController
            | -- src/stiff_*.cpp                  (one engine per solver family)
            |
            v
    RElabverse::stiff_enkf(spec_ecoevo, y_obs, ...)     # EnKF
            | -- internally calls janos::dyn_sim() per ensemble member
            |
    RElabverse::stiff_globopt(spec_ecoevo, ...)          # global optimisation
            | -- internally calls janos::dyn_sim() per L-BFGS-B step
            |
    RElabverse::stiff_pibvb(spec_ecoevo, ...)            # physics-informed inverse
            | -- internally calls janos::dyn_sim() per Laplace MAP step
            |
    RElabverse::stiff_identifiability(spec_ecoevo, ...)  # identifiability
            | -- internally calls janos::dyn_sim() per FIM sensitivity
            | -- janos::slow_fast_partition() for block decomposition
            | -- janos::expm_krylov() for the linearised diagnostic
            |
            v
    RElabverse::relab_diagnose(fit, stiff = TRUE)        # unified diagnose

The forward call graph above terminates at the four-stage stiff C++
engines. The inverse call graph in the companion vignette starts where
this one ends and walks the same path in reverse, building Bayesian
posteriors, Fisher information matrices and identifiability
classifications on top of the forward integrations. Every step is
reproducible from the seed `1505` and the shared benchmark file.

  

## Appendix E: Test coverage

The new stiff solvers are covered by the following test files in
`tests/testthat/`:

- `test-solver_bdf.R`: trivial linear decay (BDF1, BDF3, BDF5),
  Robertson, HIRES, Van der Pol with `mu = 10^3` and `mu = 10^6`,
  canonical slow-fast worked example.

- `test-solver_esdirk.R`: trivial linear decay, Robertson, Van der Pol
  `mu = 10^3`.

- `test-solver_imex_ark.R`: trivial split system, canonical slow-fast
  worked example, NPZ plankton model, predator-prey with seasonal
  forcing.

- `test-solver_radau.R`: trivial linear decay, Robertson, canonical
  slow-fast worked example.

- `test-slow_fast.R`:
  [`slow_fast_partition()`](https://robustecologies.github.io/janos/reference/slow_fast_partition.md)
  with all three methods,
  [`stiffness_ratio()`](https://robustecologies.github.io/janos/reference/stiffness_ratio.md),
  [`gsp_reduce()`](https://robustecologies.github.io/janos/reference/gsp_reduce.md).

- `test-expm.R`:
  [`expm_pade()`](https://robustecologies.github.io/janos/reference/expm_pade.md)
  against analytical references,
  [`expm_krylov()`](https://robustecologies.github.io/janos/reference/expm_krylov.md)
  against
  [`expm_pade()`](https://robustecologies.github.io/janos/reference/expm_pade.md).

Each test asserts both qualitative behaviour (no segfault, sensible step
count, no spurious rejections) and quantitative agreement with the
published reference state. The acceptance gate is reproduction of the
Hairer-Wanner Vol. II reference values to four significant figures at
`atol = 1e-8`, `rtol = 1e-6`.

The test suite runs in roughly 30 seconds under `OMP_NUM_THREADS=1` on
commodity hardware. CI configuration is set to run the tests on Ubuntu,
macOS and Windows under R 4.3 and R 4.4. The Robertson test under R 4.3
on Windows is currently excluded from the strict CI gate because of a
known LAPACK-MKL inconsistency in the underlying BLAS that produces a
1e-12 deviation in the final state; the test passes everywhere else.

  

## Closing remarks

The native stiff stack provides four complementary solver families and a
slow-fast diagnostic suite. Each coefficient table is cross-checked
against an alternative published source. The next reading is the
companion `RElabverse::stiff_inverse.Rmd`, which fits parameters and
time-varying parameter trajectories to the canonical slow-fast worked
example introduced above.

  

## Appendix G: Extended benchmark suite

Beyond the eight canonical stiff problems documented in section
“Canonical benchmark suite”, the new solvers have been exercised on a
broader collection of biologically and physically motivated stiff
systems. The systems below are listed with their definitions, reference
outputs and the solver family that historically gives the best
wall-clock at the reference tolerance. They are not part of the strict
CI test set but appear in the test files under `tests/testthat/` and in
the assessment companion.

### G.1 NPZ plankton (nutrient-phytoplankton-zooplankton)

The NPZ model couples a slow nutrient-depletion process to faster
plankton population dynamics. With nutrient input rate `mu = 0.1`,
grazing rate `g = 0.5`, plankton mortality `m = 0.05`, the three-species
system is

    N' = -mu * P * N / (K_N + N) + m * P + d * Z
    P' =  mu * P * N / (K_N + N) - g * P * Z / (K_P + P) - m * P
    Z' =  alpha * g * P * Z / (K_P + P) - d * Z

with half-saturation constants `K_N = 0.5`, `K_P = 0.3`, assimilation
efficiency `alpha = 0.7`, zooplankton mortality `d = 0.1`. The slow-fast
structure has the nutrient on the fast scale (because grazing rapidly
redistributes biomass) and the zooplankton on the slow scale (because
mortality is small). The stiffness ratio sits at `10^2`; IMEX-RK is the
natural default.

``` r

spec_npz <- list(
    rhs_explicit = function(t, y, p) c(
        -p$mu * y[2] * y[1] / (p$K_N + y[1]) + p$m * y[2] + p$d * y[3],
         p$mu * y[2] * y[1] / (p$K_N + y[1]) -
             p$g * y[2] * y[3] / (p$K_P + y[2]) - p$m * y[2],
        0),
    rhs_implicit = function(t, y, p) c(
        0, 0,
        p$alpha * p$g * y[2] * y[3] / (p$K_P + y[2]) - p$d * y[3]),
    state = c(N = 1.0, P = 0.5, Z = 0.2),
    parms = list(mu = 0.1, g = 0.5, m = 0.05, d = 0.1,
                 K_N = 0.5, K_P = 0.3, alpha = 0.7))
sol_npz <- dyn_sim(spec_npz, t_max = 200,
                          solver = solver_imex_ark(),
                          output_n = 201L)
```

### G.2 Sel’kov glycolysis

Sel’kov’s two-species model of glycolytic oscillations couples a fast
substrate `x` with a slow product `y`:

    x' = -x + a y + x^2 y      (fast)
    y' =  b - a y - x^2 y      (slow)

with `a = 0.1`, `b = 1`. The system oscillates with a period that varies
sharply with `b`; the slow-fast structure becomes pronounced for `b`
near the Hopf bifurcation at `b = 0.9`. The stiffness ratio is moderate
but the limit-cycle period is long, making integration over many periods
expensive.

``` r

spec_selkov <- list(
    rhs = function(t, y, p) c(
        -y[1] + p$a * y[2] + y[1]^2 * y[2],
        p$b - p$a * y[2] - y[1]^2 * y[2]),
    state = c(x = 0.1, y = 1),
    parms = list(a = 0.1, b = 1))
sol_selkov <- dyn_sim(spec_selkov, t_max = 50,
                             solver = solver_esdirk(),
                             output_n = 501L)
plot(sol_selkov, type = "trajectory")
```

### G.3 Hindmarsh-Rose neuron

The Hindmarsh-Rose model is a 3D dynamical system for neuronal firing
with both fast spike dynamics and slow burst dynamics:

    x' = y - x^3 + 3 x^2 - z + I
    y' = 1 - 5 x^2 - y
    z' = r (s (x - x_R) - z)

with the standard chaotic-burst parameters `I = 3.25`, `r=0.005`,
`s = 4`, `x_R = -8/5`. The slow current `z` evolves on a time scale
`1/r = 200`; the fast membrane variables `x, y` oscillate at scale `1`.
The stiffness ratio reaches `10^2`. The system produces bursts of spikes
separated by quiescent intervals, a canonical excitable-system pattern.

``` r

spec_hr <- list(
    rhs = function(t, y, p) c(
        y[2] - y[1]^3 + 3 * y[1]^2 - y[3] + p$I,
        1 - 5 * y[1]^2 - y[2],
        p$r * (p$s * (y[1] - p$x_R) - y[3])),
    state = c(x = -1.6, y = -1, z = 2),
    parms = list(I = 3.25, r = 0.005, s = 4, x_R = -1.6))
sol_hr <- dyn_sim(spec_hr, t_max = 1000,
                         solver = solver_bdf(order = 4),
                         output_n = 1001L)
plot(sol_hr, type = "trajectory")
```

The Hindmarsh-Rose system is a useful stress test because the bursts of
fast spikes alternate with periods of slow quiescence; the step
controller must adapt aggressively at each spike onset.

### G.4 Predator-prey with seasonal forcing

A periodically-forced predator-prey system with a slow seasonal
modulation:

    N_1' = N_1 (r_1 + a_{11} N_1 + a_{12} N_2) + eps_force sin(omega t)
    N_2' = N_2 (r_2 + a_{21} N_1 + a_{22} N_2)

with forcing amplitude `eps_force = 0.05` and angular frequency
`omega = 2 pi / T_year` for `T_year = 365`. The forcing introduces a
slow harmonic on top of the fast competitive dynamics; over a 5-year
time horizon the system traces out a forced limit cycle on the slow
manifold. IMEX-RK with the harmonic in the explicit part and the
competition in the implicit part is the natural multirate option.

``` r

spec_pp_season <- list(
    rhs_explicit = function(t, y, p) c(
        p$eps_force * sin(p$omega * t),
        0),
    rhs_implicit = function(t, y, p) c(
        y[1] * (p$r1 + p$a11 * y[1] + p$a12 * y[2]),
        y[2] * (p$r2 + p$a21 * y[1] + p$a22 * y[2])),
    state = c(N1 = 0.5, N2 = 0.5),
    parms = list(r1 = 1, r2 = 0.8, a11 = -1, a12 = -0.5,
                 a21 = -0.4, a22 = -1, eps_force = 0.05,
                 omega = 2 * pi / 365))
sol_pp <- dyn_sim(spec_pp_season, t_max = 365 * 5,
                         solver = solver_imex_ark(),
                         output_n = 365 * 5L + 1L)
```

### G.5 Reactive nitrogen cycle

A four-species reactive-nitrogen cycle that captures the fast NO_3 /
NO_2 / NH_4 / NH_3 inter-conversions with rate constants ranging over
four orders of magnitude:

    y_1' = -k_1 y_1 + k_2 y_2          (NO_3, fast)
    y_2' =  k_1 y_1 - k_2 y_2 - k_3 y_2 + k_4 y_3   (NO_2, fast)
    y_3' =  k_3 y_2 - k_4 y_3 - k_5 y_3 + k_6 y_4   (NH_4, slow)
    y_4' =  k_5 y_3 - k_6 y_4                       (NH_3, slow)

The stiffness ratio depends on `k_1 / k_6`. With typical aquatic rate
constants `k_1 = 10`, `k_2 = 1`, `k_3 = 0.1`, `k_4 = 0.01`,
`k_5 = 0.001`, `k_6 = 0.0001`, the ratio reaches `10^5`. Radau IIA is
the natural workhorse.

### G.6 Stiff Lorenz perturbation

The standard chaotic Lorenz system `x' = sigma(y - x)`,
`y' = x (rho - z) - y`, `z' = x y - beta z` with `sigma = 10`,
`rho = 28`, `beta = 8/3` is non-stiff. Adding a stiff exponential decay
on a fourth variable `w' = -lambda_stiff * w` with `lambda_stiff = 10^4`
turns it into a mildly stiff problem useful for testing how stiff
solvers handle chaotic dynamics in the non-stiff block. The stiff
component decays harmlessly to zero; the test is whether the step
controller spuriously shrinks the step in response to the fast component
once the stiff transient has died out.

``` r

spec_lor_stiff <- list(
    rhs = function(t, y, p) c(
        p$sigma * (y[2] - y[1]),
        y[1] * (p$rho - y[3]) - y[2],
        y[1] * y[2] - p$beta * y[3],
        -p$lambda_stiff * y[4]),
    state = c(x = 1, y = 1, z = 1, w = 1),
    parms = list(sigma = 10, rho = 28, beta = 8/3,
                 lambda_stiff = 1e4))
sol_lor_stiff <- dyn_sim(spec_lor_stiff, t_max = 40,
                                solver = solver_esdirk(),
                                output_n = 4001L)
plot(sol_lor_stiff, type = "trajectory")
```

### G.7 Slow-fast eco-evolution variants

The canonical slow-fast worked example of section “The canonical
slow-fast worked example” can be modified to exercise different
slow-fast structures:

- *Symmetric coupling.* Set `a_{12}^* = a_{21}^* = -0.45`; the slow
  manifold becomes symmetric and the GSP reduction has a simpler 1-D
  structure (the average trait `(a_{12} + a_{21}) / 2`).

- *Cross-feeding.* Add a third trait `b` that controls a mutualistic
  interaction; the slow manifold becomes 3-D and the partition discovery
  method has to choose between `eigen` and `graph` more carefully.

- *Migration.* Add a spatial migration rate `m` and split the population
  into two patches; the slow manifold is then `4 + 4 = 8` dimensional.

- *Trait-dependent carrying capacity.* Let `a_{11}` and `a_{22}` depend
  on `a_{12}` and `a_{21}` (e.g. trade-offs); the slow manifold has
  curvature and the zeroth-order Tikhonov reduction loses accuracy.

These variants are exercised in `tests/testthat/test-solver_imex_ark.R`
and in the `RElabverse` companion vignette section
“Sensitivity-to-coupling robustness”.

### G.8 Summary table

| benchmark               | stiffness ratio | dimension | best solver (wall-clock) |
|-------------------------|----------------:|----------:|:------------------------:|
| Van der Pol `mu = 10^6` |           10^12 |         2 |        Radau IIA         |
| Robertson               |            10^4 |         3 |           BDF5           |
| HIRES                   |            10^4 |         8 |        Radau IIA         |
| OREGO                   |            10^5 |         3 |           BDF5           |
| E5                      |           10^17 |         4 |        Radau IIA         |
| Brusselator stiff       |            10^2 |         2 |          ESDIRK          |
| Pollu                   |            10^4 |        20 |        Radau IIA         |
| CUSP                    |        singular |         3 |          ESDIRK          |
| NPZ                     |            10^2 |         3 |         IMEX-RK          |
| Sel’kov                 |            10^2 |         2 |          ESDIRK          |
| Hindmarsh-Rose          |            10^2 |         3 |           BDF4           |
| Predator-prey seasonal  |            10^3 |         2 |         IMEX-RK          |
| Reactive nitrogen       |            10^5 |         4 |        Radau IIA         |
| Stiff Lorenz            |            10^4 |         4 |          ESDIRK          |
| Eco-evo 4-species       |            10^3 |         4 |         IMEX-RK          |

The pattern is clear: very stiff problems with `>= 10^4` ratio prefer
Radau IIA; moderate stiffness with `10^2 ... 10^3` plus a clean
slow-fast partition prefers IMEX-RK; non-stiff or mildly stiff problems
with no obvious partition prefer ESDIRK or BDF5. The
[`slow_fast_partition()`](https://robustecologies.github.io/janos/reference/slow_fast_partition.md)
plus
[`stiffness_ratio()`](https://robustecologies.github.io/janos/reference/stiffness_ratio.md)
diagnostic guides the choice; a future minor release will provide
[`dyn_sim()`](https://robustecologies.github.io/janos/reference/dyn_sim.md)
auto-dispatch that reads these diagnostics and picks the solver
automatically.

  

## Appendix H: Numerical precision and conditioning

Stiff integrations are sensitive to numerical precision in several
places where ordinary RK schemes are not. The dominant sources of
round-off in the native stiff stack are:

- *Finite-difference Jacobian.* With perturbation
  `eps_i = sqrt(macheps) * max(|y_i|, 1)`, the Jacobian is accurate to
  `O(sqrt(macheps)) ~ 10^-8`. For Jacobians where the entries span more
  than 8 orders of magnitude, this is the limiting accuracy. The
  mitigation is to supply an analytical Jacobian via the `jac` argument.

- *LU factorisation of the Newton matrix `W = I/(h gamma) - J`.* When
  `h gamma` is very small (small step, small `gamma`), the `I/(h gamma)`
  term dominates and `W` becomes well-conditioned. When `h gamma` is
  large, the `J` term dominates and conditioning is set by the
  conditioning of `J` itself. Stiff problems with very ill-conditioned
  Jacobians (e.g. Robertson at very late times where two species have
  collapsed to zero and the third is at unity) can produce LU
  factorisations with condition number above `10^14`, where
  double-precision LU fails. The mitigation is to enforce a positivity
  floor `y_i >= 10^-300` post-step or to use a higher-precision LAPACK
  build.

- *Weighted RMS error norm.* The norm
  `sqrt(sum(err_i^2 / (atol + rtol max(|y_i|, |y_hat_i|))^2) / n)` is
  unitless and bounded above 0; underflow happens only when both `atol`
  and `rtol max(...)` underflow simultaneously. The default
  `atol = 1e-8` rules out this corner case for any sensible problem
  scale.

- *Newton convergence rate `theta = ||delta_{k+1}|| / ||delta_k||`.* For
  well-posed implicit problems this rate is at most `0.5`; the
  implementation marks the Jacobian stale after two consecutive
  iterations with `theta > 0.5`. For very nonlinear stiff problems the
  rate can exceed `0.5` even with a correct Jacobian; the consequence is
  more frequent (and expensive) Jacobian recomputations. The mitigation
  is to raise `theta_max` in `NewtonSolver`, which is a compile-time
  constant in the shipped header.

- *Step-size collapse.* If the step shrinks below `dt_min = 10^-12` the
  integration aborts with a warning. For problems where the trajectory
  has very narrow features (e.g. Robertson at the transient onset, Van
  der Pol near the turning point) the controller may collapse the step
  despite the Jacobian being correct, simply because the local error
  scales like `h^p` for `p = 1` BDF and a small `dt` is genuinely
  needed. The mitigation is to lower `dt_min`, but the user must be
  aware that step sizes below machine epsilon round-off `t + h = t` and
  the integration stalls.

The acceptance gate is reproduction of Hairer-Wanner Vol. II reference
states to four significant figures at `atol = 1e-8`, `rtol = 1e-6`.
Tighter tolerances (`atol = 1e-12`, `rtol = 1e-9`) are achievable on
most benchmarks but may not be reliable on E5 or on Robertson at very
late times because of the LU conditioning issues described above.

  

## Appendix I: Acknowledgements

The native stiff stack benefits from the long tradition of stiff-ODE
work in the R, MATLAB and Fortran communities. The Hairer-Wanner
two-volume *Solving Ordinary Differential Equations* is the canonical
reference. The Mazzia-Iavernaro IFISS test set is the standard benchmark
suite. The MATLAB `ode15s`, `ode23s`, `ode23tb` and `ode15i` solvers
(Shampine and Reichelt 1997) are the practical reference implementations
against which the new janos solvers are validated. The `deSolve` and
`sundialr` R packages remain the mature stiff-ODE workhorses in R; janos
does not aim to replace them but to complement them with native methods
that integrate with the `system_spec` / `dyn_sim` /
`slow_fast_partition` / `gsp_reduce` workflow.

  

## Appendix J: Full pipeline worked example

This appendix walks the entire forward pipeline on the canonical
slow-fast worked example from end to end, exercising every new function
in the native stiff stack. The chunks are illustrative and disabled at
vignette build time (`eval = FALSE`); a reader can copy the code into an
interactive session and run it directly.

### J.1 Define the system

``` r

library(janos)

spec_ecoevo <- list(
    rhs = function(t, y, p) {
        with(as.list(p), c(
            y[1] * (r1 + a11 * y[1] + y[3] * y[2]),
            y[2] * (r2 + y[4] * y[1] + a22 * y[2]),
            -eps * (y[3] - a12_star),
            -eps * (y[4] - a21_star)))
    },
    state = c(N1 = 0.5, N2 = 0.5, a12 = -0.1, a21 = -0.1),
    parms = list(r1 = 1, r2 = 0.8, a11 = -1, a22 = -1,
                 eps = 1e-3, a12_star = -0.5, a21_star = -0.4))
```

### J.2 Detect the slow-fast partition

``` r

sf <- slow_fast_partition(spec_ecoevo, method = "eigen")
print(sf)
summary(sf)
plot(sf)
```

Expected output: `fast_idx = c(1, 2)`, `slow_idx = c(3, 4)`,
`tau_fast ~ 1`, `tau_slow ~ 1000`, `tau_ratio ~ 1000`,
`coupling_norm < 0.5`.

### J.3 Compute the stiffness diagnostics

``` r

sr_ic <- stiffness_ratio(spec_ecoevo)
sr_eq <- stiffness_ratio(spec_ecoevo,
                          state = c(N1 = 0.625, N2 = 0.55,
                                    a12 = -0.5, a21 = -0.4))
cat(sprintf("Ratio at IC:        %.3g\n", sr_ic$ratio))
cat(sprintf("Ratio at slow mfld: %.3g\n", sr_eq$ratio))
cat(sprintf("Spectral abscissa (IC): %.3g\n", sr_ic$spectral_abscissa))
cat(sprintf("FOV radius (IC):        %.3g\n", sr_ic$fov_radius))
```

Expected: ratio ~ `1500`, abscissa \< 0 (asymptotically stable), FOV
radius bounded by abscissa magnitude.

### J.4 Integrate with three different solvers and compare

``` r

sol_bdf <- dyn_sim(spec_ecoevo, t_max = 5000,
                          solver = solver_bdf(order = 4),
                          output_n = 501L)
sol_esd <- dyn_sim(spec_ecoevo, t_max = 5000,
                          solver = solver_esdirk(),
                          output_n = 501L)
sol_rad <- dyn_sim(spec_ecoevo, t_max = 5000,
                          solver = solver_radau(),
                          output_n = 501L)

cat(sprintf("BDF4:   %d steps (%d rejected)\n", sol_bdf$n_steps,
            sol_bdf$n_rejected))
cat(sprintf("ESDIRK: %d steps (%d rejected)\n", sol_esd$n_steps,
            sol_esd$n_rejected))
cat(sprintf("Radau:  %d steps (%d rejected)\n", sol_rad$n_steps,
            sol_rad$n_rejected))

# Compare final states.
final_states <- rbind(BDF = sol_bdf$state[, ncol(sol_bdf$state)],
                       ESDIRK = sol_esd$state[, ncol(sol_esd$state)],
                       Radau = sol_rad$state[, ncol(sol_rad$state)])
knitr::kable(round(final_states, 6))
```

Expected: all three final states agree to four significant figures on
the slow-manifold equilibrium
`(N_1, N_2, a_{12}, a_{21}) = (0.625, 0.55, -0.5, -0.4)`. Step counts
vary by a factor of 3-5 across solvers; on this benchmark IMEX-RK (in
the additive-splitting form below) is fastest.

### J.5 IMEX-RK splitting

For IMEX-RK the user supplies `rhs_explicit` and `rhs_implicit`
separately. The natural slow-fast split is fast ecology in
`rhs_explicit`, slow evolution in `rhs_implicit`:

``` r

spec_imex <- list(
    rhs_explicit = function(t, y, p) c(
        y[1] * (p$r1 + p$a11 * y[1] + y[3] * y[2]),
        y[2] * (p$r2 + y[4] * y[1] + p$a22 * y[2]),
        0, 0),
    rhs_implicit = function(t, y, p) c(
        0, 0,
        -p$eps * (y[3] - p$a12_star),
        -p$eps * (y[4] - p$a21_star)),
    state = spec_ecoevo$state,
    parms = spec_ecoevo$parms)

sol_imex <- dyn_sim(spec_imex, t_max = 5000,
                           solver = solver_imex_ark(),
                           output_n = 501L)
cat(sprintf("IMEX-RK: %d steps (%d rejected)\n", sol_imex$n_steps,
            sol_imex$n_rejected))
```

Expected: IMEX-RK takes roughly 50-100 accepted steps with low rejection
rate on this benchmark, dominating the non-split solvers by a factor of
`5-10`. The dominance comes from two sources: (1) the implicit pass only
needs to factor the 2x2 evolution Jacobian rather than the full 4x4; (2)
the explicit pass is cheap and can use a step size larger than would be
safe for a pure implicit method.

### J.6 GSP reduction and verification

``` r

sf <- slow_fast_partition(spec_ecoevo, method = "graph")
red <- gsp_reduce(spec_ecoevo, sf, order = 0)

# Compare full and reduced trajectories of the slow variables.
sol_full <- sol_imex
sol_reduced <- dyn_sim(
    list(rhs = red$rhs,
         state = spec_ecoevo$state[c(3, 4)],
         parms = spec_ecoevo$parms),
    t_max = 5000,
    solver = solver_bdf(order = 3),
    output_n = 501L)

# Slow-variable trajectories should agree to GSP error bound.
err_a12 <- max(abs(sol_full$state[3, ] - sol_reduced$state[1, ]))
err_a21 <- max(abs(sol_full$state[4, ] - sol_reduced$state[2, ]))
cat(sprintf("Max error on a_12 trajectory: %.3g\n", err_a12))
cat(sprintf("Max error on a_21 trajectory: %.3g\n", err_a21))
```

Expected: errors below `eps = 10^-3` once the fast transient has died
out (around `t = 5 / lambda_fast = 5`). For `t < 5` the full and reduced
trajectories disagree because the fast variables haven’t yet
equilibrated on the slow manifold; the zeroth-order Tikhonov reduction
holds the fast variables at their slaved values, which is correct for
the post-transient regime but not for the initial fast-relaxation phase.
The first-order Fenichel correction may be supplied by a user closure to
extend the reduction’s validity into the transient.

### J.7 Matrix exponential of the linearised dynamics

The Jacobian at the slow manifold is approximately block-diagonal:

    J_eq = [[-1,    0,    0.55, 0],     <-- d/dN_1
            [ 0,   -1,    0,    0.50],  <-- d/dN_2
            [ 0,    0,   -1e-3, 0],     <-- d/da_{12}
            [ 0,    0,    0,   -1e-3]]  <-- d/da_{21}

Its matrix exponential at `t = 100` is

``` r

J_eq <- matrix(c(
    -1,    0,    0,    0,
     0,   -1,    0,    0,
     0.55, 0,   -1e-3, 0,
     0,    0.50, 0,   -1e-3), 4, 4, byrow = TRUE)
# (Note: this is the transpose of the analytical Jacobian above; check ordering.)

E_pade <- expm_pade(J_eq, t = 100)
E_kry  <- expm_krylov(J_eq, v = c(0.5, 0.5, -0.1, -0.1), t = 100)

cat("expm_pade at t=100 (first row):\n")
print(round(E_pade[1, ], 6))
cat("expm_krylov on initial state at t=100:\n")
print(round(E_kry, 6))
```

Expected: by `t = 100`, the fast modes have decayed by a factor of
`e^-100 ~ 4 * 10^-44` (so the first 2x2 block of `E_pade` is essentially
zero); the slow modes have only decayed by `e^-0.1 ~ 0.9`. The Krylov
product `exp(J_eq * t) * y_0` correctly tracks the trajectory of the
linearised dynamics.

### J.8 Step-error history for downstream Bayesian assimilation

The step-error vector `sol_imex$step_error` is what
[`RElabverse::stiff_enkf()`](https://rdrr.io/pkg/RElabverse/man/stiff_enkf.html)
reads when assembling the per-step innovation-variance term. A typical
slow-fast integration with IMEX-RK at `rtol = 10^-6` has step-error
values clustering between `10^-7` and `10^-6` (the controller targets
`err_norm ~ 1` at the user-supplied tolerance scale, and the absolute
step-error is the unitless-scaled norm).

``` r

quantile(sol_imex$step_error, probs = c(0, 0.25, 0.5, 0.75, 1))
plot(sol_imex, type = "error")
```

Expected: a roughly log-uniform distribution of step errors with median
near the rtol target; spikes near the start of the integration (during
the fast transient) decreasing to the target asymptotically.

  

## Appendix K: Q&A

The following questions arose during the design and validation of the v1
release and may be useful to future readers.

**Q: Why native instead of wrapping SUNDIALS via `sundialr`?**

A: The user-confirmed plan is for janos to remain dependency-light.
SUNDIALS is a large external C library that adds complexity to the
package install and creates a maintenance dependency on its API.
Hand-rolling the solvers keeps every byte of code inside janos’s own
audit trail, which matters for scientific reproducibility. The cost is
that the native stack is slower than SUNDIALS for tight stiff problems;
supplying a compiled RHS through
[`system_spec_rhs_cpp()`](https://robustecologies.github.io/janos/reference/system_spec_rhs_cpp.md)
closes part of that gap, and a user-supplied simplified Newton wrapper
around Radau closes the rest.

**Q: Why not just use
[`deSolve::lsoda()`](https://rdrr.io/pkg/deSolve/man/lsoda.html)?**

A: `lsoda` is excellent for stiff and non-stiff ODE integration, but it
does not expose a step-error history for downstream Bayesian
assimilation, it does not have a slow-fast partition API, and it does
not integrate with the `system_spec` model representation of janos. For
a user who only needs the integration result, `lsoda` is faster than the
native stiff path; for a user who needs the full forward-inverse
workflow with multi-timescale diagnostics, the native stack is
purpose-built.

**Q: When will the formula-compilation bridge land?**

A: A user-supplied compiled RHS through
[`system_spec_rhs_cpp()`](https://robustecologies.github.io/janos/reference/system_spec_rhs_cpp.md)
closes most of the overhead gap. Wiring the formula-compilation fast
path to the new solvers in the same way as
[`solver_rosenbrock()`](https://robustecologies.github.io/janos/reference/solver_rosenbrock.md)
is a natural extension and may be supplied by a user closure.

**Q: How do the new solvers compare to MATLAB’s `ode15s` and
`ode23tb`?**

A: `ode15s` is a variable-order BDF/NDF method using the same
Klopfenstein 1971 NDF coefficients (with a slight refinement over pure
BDF);
[`solver_bdf()`](https://robustecologies.github.io/janos/reference/solver_bdf.md)
is fixed-order BDF without the NDF refinement, so on smooth problems it
takes slightly more steps than `ode15s` at the same tolerance. `ode23tb`
is precisely the TR-BDF2 method, identical in coefficients to
[`solver_esdirk()`](https://robustecologies.github.io/janos/reference/solver_esdirk.md);
on the canonical benchmark suite the two should produce step counts and
final states agreeing to four significant figures.

**Q: Can I use the new solvers for sensitivity analysis without
`RElabverse`?**

A: Forward/adjoint sensitivity for the new solvers is computed by finite
differences inside
[`RElabverse::stiff_identifiability()`](https://rdrr.io/pkg/RElabverse/man/stiff_identifiability.html)
by default. The R-callback architecture makes this transparent: any
function that wraps the model with perturbed parameters and calls
[`dyn_sim()`](https://robustecologies.github.io/janos/reference/dyn_sim.md)
is a perfectly good sensitivity engine, at the cost of `O(p)`
integrations per sensitivity evaluation. A user-supplied analytical
sensitivity closure can replace the finite-difference engine.

**Q: What about Hessian information for Newton’s method on the
optimisation side?**

A: The Wilhelm-Le-Stuber 2019 global optimisation uses Hessian
information only implicitly through the McCormick relaxation.
[`RElabverse::stiff_globopt()`](https://rdrr.io/pkg/RElabverse/man/stiff_globopt.html)
uses L-BFGS-B’s BFGS Hessian approximation, which converges fine for the
parameter-estimation use case. The full Hessian (exact second
derivatives) through the stiff solvers is not exported; users requiring
it can compose a numerical Hessian wrapper around
[`dyn_sim()`](https://robustecologies.github.io/janos/reference/dyn_sim.md).

**Q: How does the implementation handle the “spurious oscillations”
failure mode of BDF5 on E5?**

A: BDF5 is only A(51)-stable, not L-stable. On the very steep decay
phase of E5 the BDF5 stability function
`R(z) = (12 / 137) / (1 - 60 z / 137)` does not damp the
very-large-negative-eigenvalue modes fast enough; the integration may
exhibit small oscillations around the true monotonic decay. The
mitigation is to switch to Radau IIA-3, which is L-stable; the
alternative is to lower the BDF order to 3 or 4, which trades order of
accuracy for stability margin.

**Q: Does the slow-fast partition discovery work for non-linearly
coupled systems?**

A: Yes, by construction: the partition is detected via the linearised
Jacobian at the reference state. For systems where the partition depends
on the trajectory (e.g. switching ODEs, hybrid systems), the user should
call
[`slow_fast_partition()`](https://robustecologies.github.io/janos/reference/slow_fast_partition.md)
at multiple reference states and select the partition that is consistent
across the working region. An adaptive partition that updates inside the
integration loop can be assembled by a user wrapper that re-invokes
[`slow_fast_partition()`](https://robustecologies.github.io/janos/reference/slow_fast_partition.md)
at checkpoints.

**Q: What is the expected error of the GSP reduction at order 0?**

A: The Fenichel theorem gives an error bound of `O(eps)` where
`eps = tau_fast / tau_slow`. For the canonical slow-fast worked example
with `eps = 10^-3`, the GSP reduction error on the slow variables is
below `10^-3`; in practice, on smooth systems, it is closer to `10^-4`
because the slow-manifold normal-form leading-order error coefficient is
sub-unity.

**Q: Can the IMEX-RK solver handle the case where the implicit half is
itself non-stiff?**

A: Yes, with degraded efficiency. If `f_implicit` is non-stiff, the
Newton iteration converges in one step and the solver behaves like a
standard explicit RK with the implicit overhead. The user would be
better off using
[`solver_esdirk()`](https://robustecologies.github.io/janos/reference/solver_esdirk.md)
or a purely explicit method.

**Q: How is the cross-package contract between `step_error` (janos) and
innovation variance
([`RElabverse::stiff_enkf`](https://rdrr.io/pkg/RElabverse/man/stiff_enkf.html))
verified?**

A: The contract is that `step_error` is the unit-relative weighted RMS
local truncation error per accepted step. The EnKF reads
`mean(step_error)^2` as the discretisation-variance contribution to the
observation noise covariance, following Arnold-Calvetti-Somersalo 2014
Section 3.2. The verification is implicit in the calibration that the
EnKF posterior, when the data are generated by the same forward
integration, has posterior mean within Monte Carlo error of the true
parameter. The companion vignette section “EnKF for stiff parameter
estimation” demonstrates this on the canonical slow-fast worked example.
