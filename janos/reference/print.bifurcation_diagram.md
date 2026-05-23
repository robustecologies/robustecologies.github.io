# Print method for a bifurcation_diagram object

Prints the family, parameter range, observable, grid resolution, worker
count and Lyapunov-exponent summary of a `bifurcation_diagram` object.
When the computation was interrupted by a user signal, the number of
completed parameter values is flagged with a warning glyph.

## Usage

``` r
# S3 method for class 'bifurcation_diagram'
print(x, ...)
```

## Arguments

- x:

  A `bifurcation_diagram` object.

- ...:

  Unused.

## Value

The input `x`, invisibly.

## See also

[`bifurcation_diagram()`](https://robustecologies.github.io/janos/reference/bifurcation_diagram.md)
. constructor;
[`summary.bifurcation_diagram()`](https://robustecologies.github.io/janos/reference/summary.bifurcation_diagram.md)
. attractor statistics per quartile;
[`plot.bifurcation_diagram()`](https://robustecologies.github.io/janos/reference/plot.bifurcation_diagram.md)
. orbit diagram with Lyapunov overlay.

## Examples

``` r
if (FALSE) { # \dontrun{
bd <- bifurcation_diagram(logistic_map, par_name = "r",
                           par_range = c(2.8, 4),
                           observable = "x", n_par = 200)
print(bd)
} # }
```
