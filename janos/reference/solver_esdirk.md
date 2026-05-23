# ESDIRK TR-BDF2 solver for stiff ODE systems

Three-stage explicit-first-stage singly diagonally implicit Runge-Kutta
(TR-BDF2 of Hosea-Shampine 1996), L-stable of order 2 with an embedded
order-3 error estimator. Equivalent to MATLAB's `ode23tb`.

## Usage

``` r
solver_esdirk(
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

A solver_spec object with method "esdirk".

## Details

Butcher tableau with gamma = 2 - sqrt(2), d = gamma / 2, w = sqrt(2) / 4
(Hosea-Shampine equation (1.4); cross-checked against Hairer-Wanner Vol.
II Section IV.6). The first stage is explicit, the second is a
trapezoidal-rule implicit step at t + gamma \* h, and the third is a
BDF2 step at t + h that reuses the trapezoidal-stage RHS evaluation. The
embedded order-3 estimator uses weights b_hat = ((1 - w)/3, (3 w + 1)/3,
d/3); the local error vector is h \* (b - b_hat) . k.

Jacobian reuse is governed by the modified-Newton convergence rate; the
Jacobian is rebuilt when the rate exceeds 0.5 or when the step-size
ratio leaves the interval from 0.5 to 2.0.

## References

Hosea, M. E. and Shampine, L. F. (1996). Analysis and implementation of
TR-BDF2. *ACM Trans. Math. Software* 22(1), 105-128.
[doi:10.1145/229473.229474](https://doi.org/10.1145/229473.229474) .

Hairer, E. and Wanner, G. (1996). *Solving Ordinary Differential
Equations II*. Springer, Section IV.6. ISBN 978-3-540-60452-5.

## See also

[`solver_rosenbrock`](https://robustecologies.github.io/janos/reference/solver_rosenbrock.md),
[`solver_bdf`](https://robustecologies.github.io/janos/reference/solver_bdf.md),
[`solver_imex_ark`](https://robustecologies.github.io/janos/reference/solver_imex_ark.md),
[`solver_radau`](https://robustecologies.github.io/janos/reference/solver_radau.md),
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
sol <- dyn_sim(robertson, t_max = 1e5, solver = solver_esdirk(),
               discard_transient = 0)
plot(sol)
} # }
```
