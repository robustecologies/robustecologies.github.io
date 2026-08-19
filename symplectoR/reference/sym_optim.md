# Optimize an objective with a symplectic or accelerated-gradient method

The single entry point to every solver in symplectoR. Dispatches the
objective to the requested kernel, assembles the result into a rich S3
object, and orchestrates multi-start ensembles when `n_starts` exceeds
one.

## Usage

``` r
sym_optim(
  objective,
  x0 = NULL,
  method = c("rgd", "nag", "heavy_ball", "cm", "leapfrog", "slc_poly", "slc_expo",
    "wibisono", "qhd"),
  control = sym_control(),
  n_starts = 1L,
  starts = NULL,
  n_threads = 1L,
  keep_path = c("thin", "full", "none"),
  seed = NULL,
  verbose = FALSE
)
```

## Arguments

- objective:

  A `sym_objective` built by
  [`sym_objective()`](https://robustecologies.github.io/symplectoR/reference/sym_objective.md)
  or
  [`sym_inverse()`](https://robustecologies.github.io/symplectoR/reference/sym_inverse.md).

- x0:

  Numeric start vector of length `objective$dim`. When `NULL` and the
  objective has a finite box, the box midpoint is used; otherwise an
  error is raised. Ignored by the density-based `"qhd"` method, which
  starts from a spread state over the whole box.

- method:

  Solver. `"rgd"` is relativistic gradient descent, the master kernel of
  family A; `"nag"`, `"heavy_ball"` and `"cm"` are its documented
  parameter presets (Nesterov, second-order heavy ball, classical
  momentum); `"leapfrog"` is the dissipative presymplectic leapfrog with
  selectable damping schedule; `"slc_poly"` and `"slc_expo"` are the
  polynomial and exponential Bregman leapfrog compositions with momentum
  restarting and temporal looping; `"wibisono"` is the rate-matching
  three-sequence scheme at p = 2; `"qhd"` is the classical split-step
  Fourier simulation of quantum Hamiltonian descent (requires a finite
  box, dimension at most 3).

- control:

  A `sym_control` object; see
  [`sym_control()`](https://robustecologies.github.io/symplectoR/reference/sym_control.md)
  for the per-method parameters and defaults.

- n_starts:

  Number of independent starts. Above one, start points are drawn
  (seeded, serially) from the box or from a Gaussian cloud around `x0`,
  and the result is a `sym_ensemble`.

- starts:

  Optional matrix of start points (one row per start) overriding the
  automatic generation.

- n_threads:

  OpenMP threads for the ensemble driver. Above one requires a compiled
  objective
  ([`sym_compile()`](https://robustecologies.github.io/symplectoR/reference/sym_compile.md));
  R closures fall back to a serial loop with a one-time message. `0`
  uses all available cores.

- keep_path:

  Trace retention: `"thin"` records every `record_every`-th iterate,
  `"full"` every iterate, `"none"` only scalar traces.

- seed:

  Seed for start-point generation and QHD measurement sampling.

- verbose:

  Logical; print progress.

## Value

A `sym_fit` object (single start), or a `sym_ensemble` (multiple starts)
whose `best` field is the incumbent `sym_fit`.

## Details

All family A and family B kernels evaluate exactly one gradient per
iteration. The relativistic kernel bounds each position half-kick by
`1/sqrt(delta)` (per-iteration displacement at most `2/sqrt(delta)`), so
`delta` acts as a principled trust region; its saturation fraction is
reported in the fit diagnostics. Structure-preserving discretizations do
not improve the asymptotic rate over Nesterov's method; they enlarge the
stable step-size region and preserve the phase portrait, which is where
their practical advantage lives.

## References

Franca, G., Sulam, J., Robinson, D. P., & Vidal, R. (2020). Conformal
symplectic and relativistic optimization. *Advances in Neural
Information Processing Systems*, 33, 16916-16926.

Franca, G., Jordan, M. I., & Vidal, R. (2021). On dissipative symplectic
integration with applications to gradient-based optimization. *Journal
of Statistical Mechanics: Theory and Experiment*, 2021(4), 043402.
[doi:10.1088/1742-5468/abf5d4](https://doi.org/10.1088/1742-5468/abf5d4)

Duruisseaux, V., & Leok, M. (2023). Practical perspectives on symplectic
accelerated optimization. *Optimization Methods and Software*, 38(6),
1230-1268.
[doi:10.1080/10556788.2023.2214837](https://doi.org/10.1080/10556788.2023.2214837)

Wibisono, A., Wilson, A. C., & Jordan, M. I. (2016). A variational
perspective on accelerated methods in optimization. *Proceedings of the
National Academy of Sciences*, 113(47), E7351-E7358.
[doi:10.1073/pnas.1614734113](https://doi.org/10.1073/pnas.1614734113)

## See also

[`sym_objective()`](https://robustecologies.github.io/symplectoR/reference/sym_objective.md),
[`sym_control()`](https://robustecologies.github.io/symplectoR/reference/sym_control.md),
[`sym_inverse()`](https://robustecologies.github.io/symplectoR/reference/sym_inverse.md),
[`sym_sweep()`](https://robustecologies.github.io/symplectoR/reference/sym_sweep.md),
[`print.sym_fit()`](https://robustecologies.github.io/symplectoR/reference/print.sym_fit.md),
[`summary.sym_fit()`](https://robustecologies.github.io/symplectoR/reference/summary.sym_fit.md),
[`plot.sym_fit()`](https://robustecologies.github.io/symplectoR/reference/plot.sym_fit.md)

## Examples

``` r
if (FALSE) { # \dontrun{
obj <- sym_objective(sym_benchmark("rosenbrock", d = 10))
fit <- sym_optim(obj, x0 = rep(0, 10), method = "rgd",
                 control = sym_control("rgd", eps = 0.001, delta = 10, max_iter = 20000))
print(fit)
summary(fit)
plot(fit, type = "trace")
plot(fit, type = "dashboard")
} # }
```
