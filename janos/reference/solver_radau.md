# Radau IIA solver for stiff ODE systems

Two-stage fully implicit Radau IIA of order 3, A- and L-stable
(Hairer-Wanner Vol. II Section IV.5). Block-Newton iteration on the
coupled stage system.

## Usage

``` r
solver_radau(
  order = 3,
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

- order:

  Order tag, currently only 3 supported (2-stage Radau IIA).

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

  Force the forward-difference Jacobian. Default FALSE.

## Value

A solver_spec object with method "radau".

## Details

Coefficients are taken verbatim from Hairer-Wanner Vol. II Table 5.5: c
= (1/3, 1), A = ((5/12, -1/12), (3/4, 1/4)), b = (3/4, 1/4). The coupled
stage system is solved by block-Newton iteration on the (s n)
dimensional residual R(Z) = Z - h kron(A, I) F(Y); the full Newton
matrix I - h kron(A, J) is LU-factored once per Jacobian update. The
local error estimate is the difference between stage 2 and stage 1
increments scaled by the unit-relative norm, a cheap proxy for the
Hairer-Wanner IV.8.18 estimator.

## References

Hairer, E. and Wanner, G. (1996). *Solving Ordinary Differential
Equations II*. Springer, Sections IV.5 and IV.8. ISBN 978-3-540-60452-5.

## See also

[`solver_rosenbrock`](https://robustecologies.github.io/janos/reference/solver_rosenbrock.md),
[`solver_bdf`](https://robustecologies.github.io/janos/reference/solver_bdf.md),
[`solver_esdirk`](https://robustecologies.github.io/janos/reference/solver_esdirk.md),
[`solver_imex_ark`](https://robustecologies.github.io/janos/reference/solver_imex_ark.md),
[`dyn_sim`](https://robustecologies.github.io/janos/reference/dyn_sim.md).

## Examples

``` r
if (FALSE) { # \dontrun{
robertson <- system_spec(
    rhs = list(y1 ~ -0.04 * y1 + 1e4 * y2 * y3,
               y2 ~  0.04 * y1 - 1e4 * y2 * y3 - 3e7 * y2^2,
               y3 ~  3e7 * y2^2),
    state_names = c("y1", "y2", "y3"),
    init = c(y1 = 1, y2 = 0, y3 = 0))
sol <- dyn_sim(robertson, t_max = 1e5, solver = solver_radau(),
               discard_transient = 0)
plot(sol)
} # }
```
