# Poincare section of a continuous flow

Returns the sequence of points where a continuous trajectory crosses the
codimension-one hyperplane \\\Sigma = \\y : y\_\mathrm{var} = c\\\\ in a
prescribed direction. Crossings are detected by sign changes of
\\y\_\mathrm{var}(t) - c\\ on consecutive samples and refined by linear
interpolation to sub-sample precision. Useful to reduce a
three-dimensional chaotic flow to a two-dimensional return map on the
section and to inspect the geometry of the attractor.

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

  Numeric, the level \\c\\ of the section \\\Sigma = \\y :
  y\_\mathrm{var} = c\\\\. Default: the median of `var` over the
  attractor.

- direction:

  One of `"up"` (crossings with \\\dot y\_\mathrm{var} \> 0\\), `"down"`
  (\\\dot y\_\mathrm{var} \< 0\\) or `"both"`; restricts intersections
  to the corresponding transverse orientation.

## Value

An S3 object of class `poincare_section` with components `crossings`
(data frame of interpolated intersection points, one row per crossing,
columns `time` and every state variable), `var`, `value`, `direction`
and `other_vars` (the state names not used for the section plane).

## Details

Let \\\varphi(t, y_0)\\ denote the flow of the underlying ODE. A
Poincare section is defined by a smooth codimension-one surface \\\Sigma
\subset \mathbb{R}^n\\ transverse to \\\varphi\\; the return map \\P:
\Sigma \to \Sigma\\ sends each crossing to the next. In janos the
surface is a coordinate hyperplane \\\Sigma = \\y\_\mathrm{var} = c\\\\,
and `poincare_section` returns the crossings that underlie a
discrete-time analysis of the continuous flow (periodic orbits become
fixed points of \\P^k\\; strange attractors appear as fractal clouds on
\\\Sigma\\).

The crossings are detected by examining the sign sequence of \\g(t_i) =
y\_\mathrm{var}(t_i) - c\\ on the stored attractor. For each bracket
\\(i, i+1)\\ with \\g(t_i) g(t\_{i+1}) \le 0\\ that matches `direction`,
the crossing time is obtained by linear interpolation, \\\theta =
-g(t_i) / (g(t\_{i+1}) - g(t_i))\\, and the remaining state coordinates
are likewise interpolated. No root refinement beyond linear
interpolation is performed; for chaotic flows this is below the
discretisation error of the integrator.

## References

Poincare, H. (1892). *Les methodes nouvelles de la mecanique celeste*.
Gauthier-Villars, Paris.

Guckenheimer, J., & Holmes, P. (1983). *Nonlinear Oscillations,
Dynamical Systems, and Bifurcations of Vector Fields*. Springer.
[doi:10.1007/978-1-4612-1140-2](https://doi.org/10.1007/978-1-4612-1140-2)

Strogatz, S. H. (2015). *Nonlinear Dynamics and Chaos* (2nd ed.).
Westview Press. ISBN: 978-0813349107.

## See also

[`print.poincare_section()`](https://robustecologies.github.io/janos/reference/print.poincare_section.md)
— compact header;
[`summary.poincare_section()`](https://robustecologies.github.io/janos/reference/summary.poincare_section.md)
— crossing statistics;
[`plot.poincare_section()`](https://robustecologies.github.io/janos/reference/plot.poincare_section.md)
— scatter of crossings on the section;
[`lyapunov_spectrum()`](https://robustecologies.github.io/janos/reference/lyapunov_spectrum.md)
— complementary chaos diagnostic;
[`zero_one_test()`](https://robustecologies.github.io/janos/reference/zero_one_test.md)
— scalar-observable test for chaos;
[`correlation_dimension()`](https://robustecologies.github.io/janos/reference/correlation_dimension.md)
— fractal dimension of the section cloud.

## Examples

``` r
if (FALSE) { # \dontrun{
run <- dyn_sim(rossler, t_max = 500, discard_transient = 100,
               solver = solver_rk45())
ps <- poincare_section(run, var = "z", value = 0.1)
print(ps)
summary(ps)
plot(ps)
} # }
```
