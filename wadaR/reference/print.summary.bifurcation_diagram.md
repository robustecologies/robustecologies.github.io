# Print method for a summary.bifurcation_diagram object

Formats the bifurcation-diagram summary as a tabular report:
per-quartile statistics (mean branch multiplicity, observable
percentiles, bifurcation-point density) followed by the global chaotic
fraction when Lyapunov exponents were computed.

## Usage

``` r
# S3 method for class 'summary.bifurcation_diagram'
print(x, ...)
```

## Arguments

- x:

  A `summary.bifurcation_diagram` object.

- ...:

  Unused.

## Value

The input `x`, invisibly.

## See also

[`summary.bifurcation_diagram()`](https://robustecologies.github.io/janos/reference/summary.bifurcation_diagram.md)
— summary constructor;
[`bifurcation_diagram()`](https://robustecologies.github.io/janos/reference/bifurcation_diagram.md)
— upstream constructor.

## Examples

``` r
if (FALSE) { # \dontrun{
bd <- bifurcation_diagram(logistic_map, par_name = "r",
                           par_range = c(2.8, 4),
                           observable = "x", n_par = 200)
print(summary(bd))
} # }
```
