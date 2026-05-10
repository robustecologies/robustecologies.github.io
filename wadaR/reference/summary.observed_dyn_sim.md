# Summary method for an observed_dyn_sim object

Extends
[`summary.dyn_sim()`](https://robustecologies.github.io/janos/reference/summary.dyn_sim.md)
with a signal-to-noise ratio diagnostic per observed variable: the
sample standard deviation of the true trajectory, the sample standard
deviation of the observation residual and the resulting SNR in decibels.

## Usage

``` r
# S3 method for class 'observed_dyn_sim'
summary(object, ...)
```

## Arguments

- object:

  An `observed_dyn_sim` object.

- ...:

  Unused.

## Value

A list of class `summary.observed_dyn_sim`.

## See also

[`observe()`](https://robustecologies.github.io/janos/reference/observe.md),
[`print.observed_dyn_sim()`](https://robustecologies.github.io/janos/reference/print.observed_dyn_sim.md),
[`plot.observed_dyn_sim()`](https://robustecologies.github.io/janos/reference/plot.observed_dyn_sim.md).

## Examples

``` r
if (FALSE) { # \dontrun{
run <- dyn_sim(lorenz, t_max = 10, solver = solver_rk45())
obs <- observe(run, obs_model = list(type = "gaussian", sd = 0.1))
summary(obs)
} # }
```
