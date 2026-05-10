# Summary method for a bifurcation_diagram object

Reports attractor statistics of a `bifurcation_diagram` object,
tabulated across four parameter quartiles of the scan: the mean number
of attractor points per parameter value, the 5th and 95th percentiles of
the observable, the fraction of parameter values with positive leading
Lyapunov exponent (when available), and the estimated bifurcation-point
density from changes in the branch multiplicity.

## Usage

``` r
# S3 method for class 'bifurcation_diagram'
summary(object, ...)
```

## Arguments

- object:

  A `bifurcation_diagram` object.

- ...:

  Unused.

## Value

A list of class `summary.bifurcation_diagram`, invisibly.

## Details

The branch multiplicity at a parameter value is the number of distinct
attractor points recorded after transient removal, rounded to a
tolerance of `1e-3` times the observable range. The bifurcation-point
density per quartile is the fraction of parameter values where the
branch multiplicity differs from its neighbour. It is a heuristic
estimate of local complexity and should not be interpreted as a count of
codimension-one bifurcations.

## References

Strogatz, S. H. (2015). *Nonlinear Dynamics and Chaos* (2nd ed.).
Westview Press. ISBN: 978-0813349107.

Kuznetsov, Y. A. (2004). *Elements of Applied Bifurcation Theory* (3rd
ed.). Springer.
[doi:10.1007/978-1-4757-3978-7](https://doi.org/10.1007/978-1-4757-3978-7)

## See also

[`bifurcation_diagram()`](https://robustecologies.github.io/janos/reference/bifurcation_diagram.md)
— constructor;
[`print.bifurcation_diagram()`](https://robustecologies.github.io/janos/reference/print.bifurcation_diagram.md)
— compact header;
[`plot.bifurcation_diagram()`](https://robustecologies.github.io/janos/reference/plot.bifurcation_diagram.md)
— orbit diagram with Lyapunov overlay;
[`bifurcation_sweep()`](https://robustecologies.github.io/janos/reference/bifurcation_sweep.md)
— multi-parameter bifurcation sweep;
[`continuation()`](https://robustecologies.github.io/janos/reference/continuation.md)
— pseudo-arclength continuation.

## Examples

``` r
if (FALSE) { # \dontrun{
bd <- bifurcation_diagram(logistic_map, par_name = "r",
                           par_range = c(2.8, 4),
                           observable = "x", n_par = 200)
summary(bd)
} # }
```
