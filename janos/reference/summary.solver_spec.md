# Summary method for a solver_spec object

Reports the formal classification of a `solver_spec` object: the method
family (deterministic, stochastic SDE, stochastic SSA, stochastic
tau-leap, stochastic hybrid, implicit, piecewise deterministic, spatial
stochastic), the order of consistency (weak and strong where applicable
for stochastic solvers), the step-size control regime (fixed, adaptive,
event-driven), the stability characterisation in one short sentence, and
the principal control parameters of the method.

## Usage

``` r
# S3 method for class 'solver_spec'
summary(object, ...)
```

## Arguments

- object:

  A `solver_spec` object.

- ...:

  Unused.

## Value

A list of class `summary.solver_spec`, invisibly.

## Details

The stability classification is a short textual tag, not a numerical
stability region. It states whether the method is explicit or implicit,
whether its stability region is bounded or A-stable, and whether the
method is event-driven (stochastic simulation algorithms, Lewis-Shedler
thinning), which removes the notion of step size entirely. The step-size
control regime names the embedded error estimator (Dormand-Prince
embedded pair for RK4(5), Rodas3 internal error for Rosenbrock) when the
method is adaptive.

## See also

[`solver_rk45()`](https://robustecologies.github.io/janos/reference/solver_rk45.md)
and siblings . constructors;
[`print.solver_spec()`](https://robustecologies.github.io/janos/reference/print.solver_spec.md)
. compact header;
[`plot.solver_spec()`](https://robustecologies.github.io/janos/reference/plot.solver_spec.md)
. absolute-stability region.

## Examples

``` r
if (FALSE) { # \dontrun{
summary(solver_rk45())
summary(solver_rosenbrock())
} # }
```
