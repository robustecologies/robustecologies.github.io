# Concise textual report for a simulated dynamical system

Prints a one-screen summary of a simulation produced by
[`dyn_sim`](https://robustecologies.github.io/janos/reference/dyn_sim.md),
covering the model identity, solver used, integration diagnostics
(elapsed time, step counts, rejection counts or reaction counts when
applicable), state of the attractor window, and interruption status.
Intended for interactive inspection; for numerical summaries of the
attractor use
[`summary.dyn_sim`](https://robustecologies.github.io/janos/reference/summary.dyn_sim.md)
and for visualisation use
[`plot.dyn_sim`](https://robustecologies.github.io/janos/reference/plot.dyn_sim.md).

## Usage

``` r
# S3 method for class 'dyn_sim'
print(x, ...)
```

## Arguments

- x:

  A `dyn_sim` object produced by
  [`dyn_sim`](https://robustecologies.github.io/janos/reference/dyn_sim.md).

- ...:

  Additional arguments (ignored; accepted for S3 consistency).

## Value

The input `x`, returned invisibly.

## Details

Print a dyn_sim object

The output is grouped into blocks (model, integration, attractor,
diagnostics). For stochastic solvers the seed and, when available, the
total number of reactions or rejected steps are reported. For PDMP
simulations the regime counts and switching statistics are included. The
function prints via [`cat`](https://rdrr.io/r/base/cat.html) and
[`message`](https://rdrr.io/r/base/message.html) and then returns its
input invisibly to allow chaining.

## References

Chambers, J. M. (2016). *Extending R*. Chapman and Hall/CRC. ISBN
9781498775717.

## See also

[`dyn_sim`](https://robustecologies.github.io/janos/reference/dyn_sim.md)
for the constructor;
[`summary.dyn_sim`](https://robustecologies.github.io/janos/reference/summary.dyn_sim.md)
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
run <- dyn_sim(m, t_max = 20, solver = solver_rk4(),
               discard_transient = 0, verbose = FALSE)
print(run)
} # }
```
