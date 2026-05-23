# Plot method for a bifurcation_diagram object

Renders the orbit diagram: for each parameter value, the attractor
sample is scattered along the observable axis. Map systems plot direct
iterates; flows plot local maxima of the observable. When Lyapunov
exponents have been computed, a second panel overlays the leading
exponent along the parameter axis with a reference line at zero to flag
chaotic regimes.

## Usage

``` r
# S3 method for class 'bifurcation_diagram'
plot(x, title = NULL, ...)
```

## Arguments

- x:

  A `bifurcation_diagram` object.

- title:

  Optional plot title (overrides the default).

- ...:

  Unused, kept for S3 compatibility.

## Value

A ggplot or patchwork object.

## See also

[`bifurcation_diagram()`](https://robustecologies.github.io/janos/reference/bifurcation_diagram.md)
. constructor;
[`print.bifurcation_diagram()`](https://robustecologies.github.io/janos/reference/print.bifurcation_diagram.md)
. compact header;
[`summary.bifurcation_diagram()`](https://robustecologies.github.io/janos/reference/summary.bifurcation_diagram.md)
. attractor statistics per quartile;
[`bifurcation_sweep()`](https://robustecologies.github.io/janos/reference/bifurcation_sweep.md)
. multi-parameter bifurcation sweep;
[`continuation()`](https://robustecologies.github.io/janos/reference/continuation.md)
. pseudo-arclength continuation.

## Examples

``` r
if (FALSE) { # \dontrun{
bd <- bifurcation_diagram(logistic_map, par_name = "r",
                           par_range = c(2.8, 4),
                           observable = "x", n_par = 200)
plot(bd)
} # }
```
