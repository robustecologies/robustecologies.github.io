# Milstein SDE solver

Creates a solver specification for the Milstein method for Ito SDEs. The
Milstein correction adds a term 0.5 \* g(X) \* g'(X) \* (dW^2 - dt) to
the Euler-Maruyama scheme, achieving strong order 1.0 convergence
(versus 0.5 for EM). The derivative g'(X) is approximated by forward
finite difference.

## Usage

``` r
solver_milstein(dt = 0.01, seed = 42L, subsample = 1L, dg_eps = 1e-06)
```

## Arguments

- dt:

  Time step (default: 0.01).

- seed:

  Random seed for reproducible simulations (default: 42).

- subsample:

  Store every nth step (default: 1).

- dg_eps:

  Relative perturbation for finite-difference approximation of g'(X)
  (default: 1e-6). The actual perturbation for state i is dg_eps \*
  max(\|X_i\|, 1.0).

## Value

A list of class `solver_spec` with `method = "milstein"`.

## References

Milstein, G. N. (1975). Approximate integration of stochastic
differential equations. *Theory of Probability and Its Applications*,
19(3), 557-562. [doi:10.1137/1119062](https://doi.org/10.1137/1119062)

## See also

[`solver_euler_maruyama`](https://robustecologies.github.io/janos/reference/solver_euler_maruyama.md),
[`model_spec`](https://robustecologies.github.io/janos/reference/model_spec.md),
[`dyn_sim`](https://robustecologies.github.io/janos/reference/dyn_sim.md)

## Examples

``` r
if (FALSE) { # \dontrun{
gbm <- model_spec(
    rhs = list(x ~ mu * x),
    diffusion = list(x ~ sigma * x),
    state_names = "x",
    parms = list(mu = 0.05, sigma = 0.2),
    init = c(x = 1.0)
)
result <- dyn_sim(gbm, t_max = 10, solver = solver_milstein(dt = 0.01))
plot(result)
} # }
```
