# Poincare section

Returns the points where a continuous trajectory crosses the surface
\\\\y\_\mathrm{var} = \mathrm{value}\\\\ in a given direction, computed
by linear interpolation between bracketing samples. Useful to reduce a
three-dimensional chaotic flow to a two-dimensional return map on the
section.

## Usage

``` r
poincare_section(x, var, value = NULL, direction = c("up", "down", "both"))
```

## Arguments

- x:

  A `dyn_sim` object.

- var:

  Name of the state variable defining the section plane.

- value:

  Numeric, the level \\c\\ of the section \\\\y\_\mathrm{var} = c\\\\.
  Default: median of `var` over the attractor.

- direction:

  One of `"up"`, `"down"` or `"both"`; restricts intersections to the
  corresponding crossing direction.

## Value

An S3 object of class `poincare_section` with components `crossings`
(data frame of interpolated intersection points), `var`, `value`,
`direction`, `other_vars`.

## Details

Poincare section of a 3D flow

## See also

[`lyapunov_spectrum`](https://robustecologies.github.io/janos/reference/lyapunov_spectrum.md),
[`zero_one_test`](https://robustecologies.github.io/janos/reference/zero_one_test.md)

## Examples

``` r
if (FALSE) { # \dontrun{
run <- dyn_sim(rossler, t_max = 500, discard_transient = 100,
               solver = solver_rk45())
ps <- poincare_section(run, var = "z", value = 0.1)
plot(ps)
} # }
```
