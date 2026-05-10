# Summary method for an ensemble_sim object

Produces a per-state statistical summary of the terminal states across
the ensemble: mean, standard deviation, median, first and third
quartiles, and minimum and maximum. Extinction counts are also surfaced
when the simulation recorded them.

## Usage

``` r
# S3 method for class 'ensemble_sim'
summary(object, ...)
```

## Arguments

- object:

  An `ensemble_sim` object.

- ...:

  Unused.

## Value

A list of class `summary.ensemble_sim`.

## See also

[`ensemble_sim()`](https://robustecologies.github.io/janos/reference/ensemble_sim.md),
[`print.ensemble_sim()`](https://robustecologies.github.io/janos/reference/print.ensemble_sim.md),
[`plot.ensemble_sim()`](https://robustecologies.github.io/janos/reference/plot.ensemble_sim.md).

## Examples

``` r
if (FALSE) { # \dontrun{
m <- model_spec(rhs = list(x ~ -x), state_names = "x",
                parms = list(), init = c(x = 1))
summary(ensemble_sim(m, n_replicates = 20, t_max = 1,
                      solver = solver_rk45()))
} # }
```
