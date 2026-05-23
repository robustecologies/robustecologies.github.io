# Euler-Maruyama SDE solver

Creates a solver specification for the Euler-Maruyama method for Ito
stochastic differential equations of the form dX = f(X)dt + g(X)dW,
where dW is a standard Wiener process increment.

## Usage

``` r
solver_euler_maruyama(dt = 0.01, seed = 42L, subsample = 1L)
```

## Arguments

- dt:

  Time step (default: 0.01).

- seed:

  Random seed for reproducible simulations (default: 42).

- subsample:

  Store every nth step (default: 1).

## Value

A list of class `solver_spec` with `method = "euler_maruyama"`.

## See also

[`solver_milstein`](https://robustecologies.github.io/janos/reference/solver_milstein.md),
[`system_spec`](https://robustecologies.github.io/janos/reference/system_spec.md),
[`dyn_sim`](https://robustecologies.github.io/janos/reference/dyn_sim.md)

## Examples

``` r
if (FALSE) { # \dontrun{
gbm <- system_spec(
    rhs = list(x ~ mu * x),
    diffusion = list(x ~ sigma * x),
    state_names = "x",
    parms = list(mu = 0.05, sigma = 0.2),
    init = c(x = 1.0)
)
result <- dyn_sim(gbm, t_max = 10, solver = solver_euler_maruyama(dt = 0.01))
plot(result)
} # }
```
