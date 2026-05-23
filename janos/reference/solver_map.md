# Discrete map solver

Creates a solver specification for iterating a discrete map y(n+1) =
F(y(n), parms). The `t_max` argument to `dyn_sim` is reinterpreted as
the number of iterations, and `discard_transient` as the number of
initial iterations to discard.

## Usage

``` r
solver_map(subsample = 1L)
```

## Arguments

- subsample:

  Store every nth iteration (default: 1).

## Value

A list of class `solver_spec` with `method = "map"`.

## See also

[`system_spec`](https://robustecologies.github.io/janos/reference/system_spec.md),
[`dyn_sim`](https://robustecologies.github.io/janos/reference/dyn_sim.md)

## Examples

``` r
if (FALSE) { # \dontrun{
logistic <- system_spec(
    map = list(x ~ r * x * (1 - x)),
    state_names = "x",
    parms = list(r = 3.9),
    init = c(x = 0.1)
)
result <- dyn_sim(logistic, t_max = 10000, solver = solver_map(),
                  discard_transient = 500)
plot(result)
} # }
```
