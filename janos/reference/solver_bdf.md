# Variable-order BDF solver for stiff ODE systems

Backward differentiation formula of fixed order (1 to 5) with adaptive
step size and modified-Newton iteration. Order is selected at
construction and held constant across the integration.

## Usage

``` r
solver_bdf(
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

  BDF order; integer in {1, 2, 3, 4, 5}. Default 3.

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

  Force the forward-difference Jacobian even when the model carries an
  analytical Jacobian. Default FALSE.

## Value

A solver_spec object with method "bdf".

## Details

BDF coefficients are tabulated from Hairer-Wanner Vol. II, Section
III.1, Table 1.1. The Newton iteration uses an analytical or
forward-difference Jacobian and reuses the Jacobian across steps when
the step-size ratio remains within the interval from 0.5 to 2.0 and the
modified-Newton convergence rate stays below 0.5. Step-size control
follows the Gustafsson PI predictive scheme (safety 0.9, alpha 0.7, beta
0.4). The local error estimate compares the order-k BDF solution against
an order-(k-1) reference following Shampine 1980.

Order 5 is the upper safe limit before zero-stability fails on stiff
problems; the constructor errors out for orders below 1 or above 5. The
integration cold-starts with k-1 implicit-Euler steps to populate the
BDF(k) history.

## References

Hairer, E. and Wanner, G. (1996). *Solving Ordinary Differential
Equations II: Stiff and Differential-Algebraic Problems*. Springer. ISBN
978-3-540-60452-5.

Brown, P. N., Byrne, G. D. and Hindmarsh, A. C. (1989). VODE: A
variable-coefficient ODE solver. *SIAM J. Sci. Stat. Comput.* 10(5),
1038-1051. [doi:10.1137/0910062](https://doi.org/10.1137/0910062) .

Shampine, L. F. (1980). Implementation of implicit formulas for the
solution of ODEs. *SIAM J. Sci. Stat. Comput.* 1(1), 103-118.
[doi:10.1137/0901005](https://doi.org/10.1137/0901005) .

## See also

[`solver_rosenbrock`](https://robustecologies.github.io/janos/reference/solver_rosenbrock.md),
[`solver_esdirk`](https://robustecologies.github.io/janos/reference/solver_esdirk.md),
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
sol <- dyn_sim(robertson, t_max = 1e5, solver = solver_bdf(order = 3),
               discard_transient = 0)
plot(sol)
} # }
```
