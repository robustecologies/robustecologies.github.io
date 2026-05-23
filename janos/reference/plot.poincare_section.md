# Plot method for a poincare_section object

Renders the crossings of a `poincare_section` object as a scatter plot
on two of the non-sectioning state variables. The resulting cloud is the
image of the attractor on the section plane: periodic orbits appear as
discrete points, tori as closed curves, strange attractors as fractal
clouds.

## Usage

``` r
# S3 method for class 'poincare_section'
plot(x, vars = NULL, title = NULL, ...)
```

## Arguments

- x:

  A `poincare_section` object.

- vars:

  Character vector of length 2 giving the two non-sectioning state
  variables to plot on the x and y axes. Defaults to the first two
  entries of `x$other_vars`.

- title:

  Optional plot title (overrides the default).

- ...:

  Unused, kept for S3 compatibility.

## Value

A ggplot object.

## See also

[`poincare_section()`](https://robustecologies.github.io/janos/reference/poincare_section.md)
. constructor;
[`print.poincare_section()`](https://robustecologies.github.io/janos/reference/print.poincare_section.md)
. compact header;
[`summary.poincare_section()`](https://robustecologies.github.io/janos/reference/summary.poincare_section.md)
. crossing statistics.

## Examples

``` r
if (FALSE) { # \dontrun{
run <- dyn_sim(rossler, t_max = 500, discard_transient = 100)
plot(poincare_section(run, var = "z", value = 0.1))
} # }
```
