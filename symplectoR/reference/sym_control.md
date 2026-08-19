# Control parameters for symplectic optimizers

Assembles and validates the control list consumed by
[`sym_optim()`](https://robustecologies.github.io/symplectoR/reference/sym_optim.md).
Any parameter not supplied falls back to the per-method default; unknown
parameters for the chosen method raise an error rather than being
silently ignored.

## Usage

``` r
sym_control(method = NULL, ...)
```

## Arguments

- method:

  Optional method name; when supplied, parameters are validated against
  that method immediately. When `NULL`, validation is deferred to
  [`sym_optim()`](https://robustecologies.github.io/symplectoR/reference/sym_optim.md),
  which knows the method.

- ...:

  Named control parameters. Shared across methods: `max_iter` (iteration
  budget, default 1000), `tol_grad` (gradient-norm stopping tolerance,
  default 1e-8), `tol_f` (objective-change stopping tolerance, default 0
  meaning disabled), `record_every` (trace thinning stride), `verbose`.
  Method-specific parameters:

  `rgd`

  :   `eps` (learning rate), `mu` (momentum factor in (0,1)), `delta`
      (relativistic factor; each of the two position half-kicks is
      bounded by `1/sqrt(delta)`, so the per-iteration displacement is
      bounded by `2/sqrt(delta)`; 0 disables), `alpha` (symplecticity
      interpolation in `[0,1]`; 1 is conformal symplectic and
      recommended).

  `nag`, `heavy_ball`, `cm`

  :   `eps`, `mu` only; these are documented parameter presets of the
      relativistic kernel.

  `leapfrog`

  :   `h` (step), `gamma` (damping rate), `damping` (schedule:
      `"constant"` is eta(t) = gamma t, `"nesterov"` is gamma log(1 +
      t), `"mixed"` adds `gamma2 * t^delta_exp`), `t0`.

  `slc_poly`

  :   `p` (polynomial order, default 6), `C`, `h`, `restart` (gradient
      momentum restart), `temporal_loop`, `loop_beta`, `loop_eps`.

  `slc_expo`

  :   `eta` (rate, default 0.01), `C`, `h`, and the same safeguards.

  `wibisono`

  :   `eps` (step), `N` (regularization factor, above 1), `C`
      (estimate-sequence constant, at most `1/(8 N)` for the p = 2
      rate).

  `qhd`

  :   `lambda` (schedule: `"cubic"`, `"expo"`, `"onethird"`), `N_grid`
      (points per dimension), `T_evol`, `K` (steps), `mu_qhd`,
      `alpha_nc`, `n_samples`, `L_box`.

## Value

An object of class `sym_control`: the validated parameter list with the
method name attached when supplied.

## Details

Family B defaults follow the tuning study of Duruisseaux and Leok
(2023): fix the family order (`p = 6` or `eta = 0.01`) and tune only
`(C, h)`, starting from `(0.1, 0.01)` for the polynomial family and
`(1, 4)` for the exponential family; the convergent region of the
`(C, h)` plane is nearly dimension independent, so tuning on a small
analogue transfers.

## References

Duruisseaux, V., & Leok, M. (2023). Practical perspectives on symplectic
accelerated optimization. *Optimization Methods and Software*, 38(6),
1230-1268.
[doi:10.1080/10556788.2023.2214837](https://doi.org/10.1080/10556788.2023.2214837)

## See also

[`sym_optim()`](https://robustecologies.github.io/symplectoR/reference/sym_optim.md),
[`sym_sweep()`](https://robustecologies.github.io/symplectoR/reference/sym_sweep.md),
[`sym_objective()`](https://robustecologies.github.io/symplectoR/reference/sym_objective.md)

## Examples

``` r
if (FALSE) { # \dontrun{
ctrl <- sym_control("rgd", eps = 0.05, mu = 0.95, delta = 4, max_iter = 5000)
fit <- sym_optim(sym_objective(sym_benchmark("rosenbrock", d = 10)),
                 x0 = rep(0, 10), method = "rgd", control = ctrl)
} # }
```
