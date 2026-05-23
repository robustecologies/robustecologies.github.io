# Print method for a summary.observed_dyn_sim object

Formats the observed-simulation summary as a compact report including
the base `dyn_sim` summary followed by the SNR table per observed
variable.

## Usage

``` r
# S3 method for class 'summary.observed_dyn_sim'
print(x, ...)
```

## Arguments

- x:

  A `summary.observed_dyn_sim` object.

- ...:

  Unused.

## Value

The input `x`, invisibly.

## See also

[`summary.observed_dyn_sim()`](https://robustecologies.github.io/janos/reference/summary.observed_dyn_sim.md),
[`observe()`](https://robustecologies.github.io/janos/reference/observe.md).

## Examples

``` r
if (FALSE) { # \dontrun{
run <- dyn_sim(lorenz, t_max = 10, solver = solver_rk45())
obs <- observe(run, obs_model = list(type = "gaussian", sd = 0.1))
print(summary(obs))
} # }
```
