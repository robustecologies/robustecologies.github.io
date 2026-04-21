# Plot method for a lyapunov_function object

Visualises a Lyapunov function \\V(x)\\ and its orbital derivative
\\\dot V(x)\\ for two-dimensional systems. Six plot types are available,
ranging from contour maps to 3D surfaces and phase portraits with
integrated trajectories.

## Usage

``` r
# S3 method for class 'lyapunov_function'
plot(
  x,
  type = c("level_sets", "landscape", "gradient_field", "phase_portrait", "V_dot",
    "comparison", "cobweb", "iterate_decay", "generator_field", "ensemble_decay",
    "noise_ellipse", "lmi_spectrum", "delay_margin", "regime_lmi", "switched_trajectory",
    "drift_grid", "fluid_vs_ctmc", "energy_decay", "gradient_field_check", "profile"),
  n_grid = 60L,
  xlim = NULL,
  ylim = NULL,
  n_trajectories = 12L,
  ...
)
```

## Arguments

- x:

  A `lyapunov_function` object.

- type:

  Character, one of `"level_sets"` (default), `"landscape"`,
  `"gradient_field"`, `"phase_portrait"`, `"V_dot"` or `"comparison"`.

- n_grid:

  Integer, number of grid points per axis. Defaults to 60.

- xlim, ylim:

  Numeric vectors of length 2 giving the plotting range. Defaults to a
  region around `x*`.

- n_trajectories:

  Integer, number of sample trajectories for `"phase_portrait"`.
  Defaults to 12.

- ...:

  Unused, kept for S3 compatibility.

## Value

A ggplot or plotly object.

## Details

`"level_sets"` shows filled contours of \\V(x)\\ with the RElab
sequential palette, the equilibrium as a red point with white halo, and
(for Goh/MacArthur) the positive orthant boundary. `"landscape"`
produces an interactive 3D plotly surface of \\V(x)\\ colored by \\\dot
V(x)\\. `"gradient_field"` overlays the vector field \\f(x)\\ as arrows
colored by \\\dot V\\ magnitude using the RElab diverging palette.
`"phase_portrait"` integrates sample trajectories and overlays them on
\\V\\ contours. `"V_dot"` plots the orbital derivative \\\dot V(x)\\
alone, with the zero-contour as a thick black boundary delineating the
effective domain of attraction. `"comparison"` places \\V\\ and \\\dot
V\\ side by side via patchwork.

## References

Lyapunov, A. M. (1892). *The general problem of the stability of
motion*. Translation in: Int. J. Control, 55(3), 531-773 (1992).
[doi:10.1080/00207179208934253](https://doi.org/10.1080/00207179208934253)

Gilpin, M. E. (1974). A Liapunov function for competition communities.
*Journal of Theoretical Biology*, 44(1), 35-48.
[doi:10.1016/S0022-5193(74)80028-7](https://doi.org/10.1016/S0022-5193%2874%2980028-7)

## See also

[`lyapunov_function()`](https://robustecologies.github.io/janos/reference/lyapunov_function.md),
[`print.lyapunov_function()`](https://robustecologies.github.io/janos/reference/print.lyapunov_function.md),
[`summary.lyapunov_function()`](https://robustecologies.github.io/janos/reference/summary.lyapunov_function.md)

## Examples

``` r
if (FALSE) { # \dontrun{
A <- matrix(c(-1, 0.5, 0, -2), 2, 2)
plot(lyapunov_function(A))
plot(lyapunov_function(A), type = "landscape")
plot(lyapunov_function(A), type = "comparison")

r <- c(1, 0.5); alpha <- matrix(c(-1, -0.3, -0.2, -1), 2, 2)
plot(lyapunov_function(list(r = r, alpha = alpha)), type = "phase_portrait")
} # }
```
