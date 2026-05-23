# Additive IMEX-RK ARS(2,3,2) for slow-fast stiff systems

Two-stage implicit, three-stage explicit additive Runge-Kutta of order 2
with an order-1 embedded estimator (Ascher-Ruuth-Spiteri 1997). The
model must supply both an explicit (`rhs_explicit`) and an implicit
(`rhs_implicit`) right-hand side. The total derivative is f(t, y) =
rhs_explicit(t, y) + rhs_implicit(t, y); the implicit half is handled by
the SDIRK pass and the explicit half by the three-stage ERK pass.

## Usage

``` r
solver_imex_ark(
  method = "ars232",
  dt_init = 0.01,
  atol = 1e-08,
  rtol = 1e-06,
  dt_min = 1e-12,
  dt_max = 1,
  max_steps = 1e+07,
  output_n = 201L,
  use_fd_jac = FALSE
)
```

## Arguments

- method:

  Variant tag, currently only `"ars232"`.

- dt_init:

  Initial step size. Default 0.01.

- atol:

  Absolute tolerance. Default 1e-8.

- rtol:

  Relative tolerance. Default 1e-6.

- dt_min:

  Lower bound on the step size. Default 1e-12.

- dt_max:

  Upper bound on the step size. Default 1.0.

- max_steps:

  Maximum total accepted steps. Default 1e7.

- output_n:

  Number of evenly spaced output points returned to R. Default 201.

- use_fd_jac:

  Force the forward-difference implicit Jacobian. Default FALSE.

## Value

A solver_spec object with method "imex_ark".

## Details

For slow-fast systems the canonical split assigns the stiff component
(large negative eigenvalues, source of step-size restriction under
explicit integration) to `rhs_implicit`, and the non-stiff component
(cheap to evaluate but otherwise unconstrained) to `rhs_explicit`. The
implicit pass works on the small stiff block while the explicit pass
tracks the cheap block at low cost per step. Coefficients are taken
verbatim from Ascher-Ruuth-Spiteri 1997 Table 1: gamma = (2 - sqrt(2))/2
and delta = 1 - 1/(2 gamma).

The implicit Jacobian, if supplied, must be the Jacobian of the implicit
half only. If absent, a forward-difference Jacobian of `rhs_implicit` is
computed.

## References

Ascher, U. M., Ruuth, S. J. and Spiteri, R. J. (1997). Implicit-explicit
Runge-Kutta methods for time-dependent partial differential equations.
*Appl. Numer. Math.* 25(2-3), 151-167.
[doi:10.1016/S0168-9274(97)00056-1](https://doi.org/10.1016/S0168-9274%2897%2900056-1)
.

Kennedy, C. A. and Carpenter, M. H. (2003). Additive Runge-Kutta schemes
for convection-diffusion-reaction equations. *Appl. Numer. Math.*
44(1-2), 139-181.
[doi:10.1016/S0168-9274(02)00138-1](https://doi.org/10.1016/S0168-9274%2802%2900138-1)
.

## See also

[`solver_bdf`](https://robustecologies.github.io/janos/reference/solver_bdf.md),
[`solver_esdirk`](https://robustecologies.github.io/janos/reference/solver_esdirk.md),
[`solver_radau`](https://robustecologies.github.io/janos/reference/solver_radau.md),
[`dyn_sim`](https://robustecologies.github.io/janos/reference/dyn_sim.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Slow-fast split: cheap non-stiff block explicit, stiff block implicit.
fE <- function(t, y, p) c(y[1] * (p$r - p$a * y[2]), 0, 0)
fI <- function(t, y, p) c(0, -p$eps * (y[2] - y[3]),
                          -p$eps * (y[3] - 0.5))
spec <- system_spec(
    rhs_explicit = fE, rhs_implicit = fI,
    state_names = c("N", "x", "x_star"),
    init = c(N = 1, x = 0.1, x_star = 0.5),
    parms = list(r = 0.5, a = 0.4, eps = 1e-3))
sol <- dyn_sim(spec, t_max = 5000, solver = solver_imex_ark(),
               discard_transient = 0)
plot(sol)
} # }
```
