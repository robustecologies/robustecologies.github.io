# Numerical summary of a simulated dynamical system

Computes and returns descriptive statistics on the attractor portion of
a simulation produced by
[`dyn_sim`](https://robustecologies.github.io/janos/reference/dyn_sim.md).
For non-spatial models the output reports, for every state variable, the
minimum, maximum, mean, standard deviation and median over the
post-transient window. For PDE and RDME models, species-level statistics
are first reduced over the spatial dimension and then over time. Useful
for quickly characterising the long-term behaviour (fixed point, limit
cycle, noise level, invariant density moments).

## Usage

``` r
# S3 method for class 'dyn_sim'
summary(object, ...)
```

## Arguments

- object:

  A `dyn_sim` object produced by
  [`dyn_sim`](https://robustecologies.github.io/janos/reference/dyn_sim.md).

- ...:

  Additional arguments (ignored; accepted for S3 consistency).

## Value

A list of class `summary.dyn_sim` with components for the per-state
numerical summaries, the names of the state variables, the solver label
and the attractor window. The object has its own `print` method.

## Details

Summary statistics for a dyn_sim object

Attractor extraction follows the convention of
[`dyn_sim`](https://robustecologies.github.io/janos/reference/dyn_sim.md):
all time points with \\t \ge\\ `discard_transient` are retained. When no
such points exist the function falls back to the full trajectory and
issues a warning.

## References

Strogatz, S. H. (2015). *Nonlinear Dynamics and Chaos*, 2nd edition.
Westview Press. ISBN 9780813349107.

## See also

[`dyn_sim`](https://robustecologies.github.io/janos/reference/dyn_sim.md)
for the constructor;
[`print.dyn_sim`](https://robustecologies.github.io/janos/reference/print.dyn_sim.md)
and
[`plot.dyn_sim`](https://robustecologies.github.io/janos/reference/plot.dyn_sim.md)
for the companion S3 methods;
[`system_spec`](https://robustecologies.github.io/janos/reference/system_spec.md)
for defining the underlying system.

## Examples

``` r
if (FALSE) { # \dontrun{
m <- system_spec(
    rhs = list(x ~ -0.1 * x + y, y ~ -x - 0.1 * y),
    state_names = c("x", "y"),
    parms = list(),
    init  = c(x = 1, y = 0)
)
run <- dyn_sim(m, t_max = 40, solver = solver_rk4(),
               discard_transient = 10, verbose = FALSE)
summary(run)
} # }
```
