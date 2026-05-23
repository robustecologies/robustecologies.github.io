# Print method for an ensemble_sim object

One-screen header for an `ensemble_sim` object: number of replicates,
solver, simulation duration, parallel backend and worker count, elapsed
time, per-state mean and standard deviation of the terminal state, and
extinction-event counts when the metadata carries them.

## Usage

``` r
# S3 method for class 'ensemble_sim'
print(x, ...)
```

## Arguments

- x:

  An `ensemble_sim` object.

- ...:

  Unused.

## Value

The input `x`, invisibly.

## See also

[`ensemble_sim()`](https://robustecologies.github.io/janos/reference/ensemble_sim.md),
[`summary.ensemble_sim()`](https://robustecologies.github.io/janos/reference/summary.ensemble_sim.md),
[`plot.ensemble_sim()`](https://robustecologies.github.io/janos/reference/plot.ensemble_sim.md).

## Examples

``` r
if (FALSE) { # \dontrun{
m <- system_spec(rhs = list(x ~ -x), state_names = "x",
                parms = list(), init = c(x = 1))
ens <- ensemble_sim(m, n_replicates = 20, t_max = 1,
                     solver = solver_rk45())
print(ens)
} # }
```
