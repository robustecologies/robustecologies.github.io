# Summary method for a poincare_section object

Reports quantitative statistics of a `poincare_section` object: the
number of crossings, the mean and median inter-crossing time (which
approximates the local period of the flow on the section), the ranges
and medians of each non-sectioning state variable on the section, and
the number of transverse crossings in each direction when
`direction = "both"`.

## Usage

``` r
# S3 method for class 'poincare_section'
summary(object, ...)
```

## Arguments

- object:

  A `poincare_section` object.

- ...:

  Unused.

## Value

A list of class `summary.poincare_section` with the above statistics,
invisibly.

## Details

The inter-crossing time \\\Delta t_k = t\_{k+1} - t_k\\ is the time
between consecutive intersections with the section plane and is an
empirical estimate of the return time of the Poincare map. For a
periodic orbit of period \\T\\ crossing the section \\m\\ times per
revolution, \\\Delta t_k \to T / m\\; for chaotic flows the distribution
of \\\Delta t_k\\ is typically multi-modal. The ranges and medians of
the other state variables are useful to locate the attractor on the
section plane before plotting.

## References

Guckenheimer, J., & Holmes, P. (1983). *Nonlinear Oscillations,
Dynamical Systems, and Bifurcations of Vector Fields*. Springer.
[doi:10.1007/978-1-4612-1140-2](https://doi.org/10.1007/978-1-4612-1140-2)

## See also

[`poincare_section()`](https://robustecologies.github.io/janos/reference/poincare_section.md)
— constructor;
[`print.poincare_section()`](https://robustecologies.github.io/janos/reference/print.poincare_section.md)
— compact header;
[`plot.poincare_section()`](https://robustecologies.github.io/janos/reference/plot.poincare_section.md)
— scatter of crossings on the section.

## Examples

``` r
if (FALSE) { # \dontrun{
run <- dyn_sim(rossler, t_max = 500, discard_transient = 100,
               solver = solver_rk45())
summary(poincare_section(run, var = "z", value = 0.1))
} # }
```
