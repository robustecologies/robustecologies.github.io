# Print method for an observed_dyn_sim object

Extends
[`print.dyn_sim()`](https://robustecologies.github.io/janos/reference/print.dyn_sim.md)
with a one-screen description of the observation layer: the noise model,
its parameters, the set of observed variables and the seed used for
reproducibility.

## Usage

``` r
# S3 method for class 'observed_dyn_sim'
print(x, ...)
```

## Arguments

- x:

  An `observed_dyn_sim` object.

- ...:

  Additional arguments passed to
  [`print.dyn_sim()`](https://robustecologies.github.io/janos/reference/print.dyn_sim.md).

## Value

The input `x`, invisibly.

## See also

[`observe()`](https://robustecologies.github.io/janos/reference/observe.md),
[`summary.observed_dyn_sim()`](https://robustecologies.github.io/janos/reference/summary.observed_dyn_sim.md),
[`plot.observed_dyn_sim()`](https://robustecologies.github.io/janos/reference/plot.observed_dyn_sim.md).

## Examples

``` r
if (FALSE) { # \dontrun{
run <- dyn_sim(lorenz, t_max = 10, solver = solver_rk45())
obs <- observe(run, obs_model = list(type = "gaussian", sd = 0.1))
print(obs)
} # }
```
