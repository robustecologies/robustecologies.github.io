# DDE solver with fixed-step RK4 and history buffer

Creates a solver specification for delay differential equations using
fixed-step fourth-order Runge-Kutta integration with Hermite cubic
interpolation of the history buffer. The history buffer stores (t, y,
dydt) triples in a circular buffer for efficient delay lookups.

## Usage

``` r
solver_dde(dt = 0.1, subsample = 1L, history_resolution = 1L)
```

## Arguments

- dt:

  Time step (default: 0.1).

- subsample:

  Store every nth step (default: 1).

- history_resolution:

  Number of extra history points per delay interval for interpolation
  accuracy (default: 1). Increasing this value improves interpolation
  quality but uses more memory.

## Value

A list of class `solver_spec` with `method = "dde"`.

## See also

[`system_spec`](https://robustecologies.github.io/janos/reference/system_spec.md),
[`dyn_sim`](https://robustecologies.github.io/janos/reference/dyn_sim.md)

## Examples

``` r
if (FALSE) { # \dontrun{
mg <- system_spec(
    rhs = list(x ~ (a * x_delay) / (1 + x_delay^c) - b * x),
    delays = list(x_delay = list(state = "x", tau = 17.0)),
    state_names = "x",
    parms = list(a = 0.2, b = 0.1, c = 10.0),
    init = c(x = 0.9)
)
result <- dyn_sim(mg, t_max = 1000, solver = solver_dde(dt = 0.1),
                  discard_transient = 500)
plot(result)
} # }
```
