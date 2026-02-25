# Compute Newton fractal basins

Computes basins of attraction for the Newton-Raphson iteration applied
to the polynomial \\z^n - 1 = 0\\. Uses a high-performance Rcpp
implementation with OpenMP parallelization.

## Usage

``` r
compute_newton_basins(
  n_roots = 3,
  x_range = c(-2, 2),
  y_range = c(-2, 2),
  resolution = 500,
  max_iter = 100,
  tolerance = 1e-06,
  n_cores = NULL
)
```

## Arguments

- n_roots:

  Integer. Degree of the polynomial (number of roots of unity). Must be
  at least 3 for Wada basin analysis. Default is 3.

- x_range:

  Numeric vector of length 2. Range of real part \\\[\text{Re}\_{\min},
  \text{Re}\_{\max}\]\\. Default is `c(-2, 2)`.

- y_range:

  Numeric vector of length 2. Range of imaginary part
  \\\[\text{Im}\_{\min}, \text{Im}\_{\max}\]\\. Default is `c(-2, 2)`.

- resolution:

  Integer. Number of grid points in each dimension. Total grid size is
  `resolution x resolution`. Default is 500.

- max_iter:

  Integer. Maximum number of Newton iterations per starting point.
  Default is 100.

- tolerance:

  Numeric. Convergence tolerance. A point is considered converged to
  root \\\zeta_k\\ when \\\|z - \zeta_k\| \< \text{tolerance}\\. Default
  is `1e-6`.

## Value

An object of class `wada_basins` (list) containing:

- basins:

  Integer matrix (`resolution x resolution`). Basin assignments where
  `basins[i,j] = k` means the point converged to the \\k\\-th root.
  Value 0 indicates non-convergence within `max_iter`.

- x_grid:

  Numeric vector of real part grid values.

- y_grid:

  Numeric vector of imaginary part grid values.

- n_attractors:

  Integer. Equal to `n_roots`.

- unclassified:

  Integer. Count of points that did not converge.

- x_range, y_range:

  The input ranges.

- resolution:

  The input resolution.

## Details

The Newton-Raphson method for finding roots of \\f(z) = z^n - 1\\
iterates: \$\$z\_{k+1} = z_k - \frac{z_k^n - 1}{n z_k^{n-1}} = z_k \cdot
\frac{(n-1)z_k^n + 1}{n z_k^n}\$\$

The polynomial \\z^n - 1 = 0\\ has exactly \\n\\ roots (the \\n\\-th
roots of unity): \$\$\zeta_k = e^{2\pi i k / n} = \cos\left(\frac{2\pi
k}{n}\right) + i\sin\left(\frac{2\pi k}{n}\right)\$\$ for \\k = 0, 1,
\ldots, n-1\\. These roots are equally spaced on the unit circle.

The basin of attraction for root \\\zeta_k\\ is: \$\$B_k = \\z_0 \in
\mathbb{C} : \lim\_{m \to \infty} N^m(z_0) = \zeta_k\\\$\$ where \\N\\
is the Newton map.

**Wada property:** For \\n \geq 3\\, the Newton fractal exhibits the
Wada property. The boundary \\\partial B_k\\ of any basin coincides with
the Julia set of the Newton map, which is the same for all basins:
\$\$\partial B_1 = \partial B_2 = \cdots = \partial B_n =
\mathcal{J}(N)\$\$

This means every boundary point is simultaneously on the boundary of all
\\n\\ basins—a remarkable topological property with implications for the
fundamental unpredictability of root-finding algorithms.

The function uses parallel computation via OpenMP to efficiently
evaluate the Newton iteration across all grid points.

## References

Curry, J. H., Garnett, L., & Sullivan, D. (1983). On the iteration of a
rational function: Computer experiments with Newton's method.
*Communications in Mathematical Physics*, 91(2), 267-277.
[doi:10.1007/BF01211162](https://doi.org/10.1007/BF01211162)

Przytycki, F. (1989). Remarks on the simple connectedness of basins of
sinks for iterations of rational maps. *Banach Center Publications*,
23(1), 229-235.
[doi:10.4064/-23-1-229-235](https://doi.org/10.4064/-23-1-229-235)

Peitgen, H.-O., & Richter, P. H. (1986). *The beauty of fractals: Images
of complex dynamical systems*. Springer-Verlag.
[doi:10.1007/978-3-642-61717-1](https://doi.org/10.1007/978-3-642-61717-1)

Daza, A., Wagemakers, A., & Sanjuan, M. A. F. (2022). Classifying basins
of attraction using the basin entropy. *Chaos, Solitons & Fractals*,
159, 112112.
[doi:10.1016/j.chaos.2022.112112](https://doi.org/10.1016/j.chaos.2022.112112)

## See also

[`newton_fractal_system`](https://robustecologies.github.io/wadaR/reference/newton_fractal_system.md)
for creating a Newton system object,
[`compute_basins`](https://robustecologies.github.io/wadaR/reference/compute_basins.md)
for general basin computation,
[`wada_grid_method`](https://robustecologies.github.io/wadaR/reference/wada_grid_method.md)
for testing the Wada property.

## Examples

``` r
if (FALSE) { # \dontrun{
# ===================================================================== #
# Example 1: Classic cubic Newton fractal (z^3 - 1 = 0)                 #
# ===================================================================== #
basins3 <- compute_newton_basins(n_roots = 3, resolution = 500)
plot(basins3, title = "Newton fractal: z^3 - 1 = 0",
     colors = c("#E41A1C", "#377EB8", "#4DAF4A"))

# The three roots are at exp(2*pi*i*k/3) for k = 0, 1, 2
# i.e., at (1, 0), (-0.5, sqrt(3)/2), (-0.5, -sqrt(3)/2)

# ===================================================================== #
# Example 2: Quintic Newton fractal (z^5 - 1 = 0)                       #
# ===================================================================== #
basins5 <- compute_newton_basins(n_roots = 5, resolution = 600)
plot(basins5, title = "Newton fractal: z^5 - 1 = 0")

# ===================================================================== #
# Example 3: Zoom into boundary region                                  #
# ===================================================================== #
basins_zoom <- compute_newton_basins(n_roots = 3,
                                     x_range = c(-0.5, 0.5),
                                     y_range = c(0.5, 1.5),
                                     resolution = 800)
plot(basins_zoom, title = "Newton fractal: zoomed boundary")

# ===================================================================== #
# Example 4: High-degree polynomial for complex fractal                 #
# ===================================================================== #
basins7 <- compute_newton_basins(n_roots = 7, resolution = 500)
plot(basins7, title = "Newton fractal: z^7 - 1 = 0")

# ===================================================================== #
# Example 5: Test for Wada property                                     #
# ===================================================================== #
basins <- compute_newton_basins(n_roots = 3, resolution = 400)
result <- wada_grid_method(basins)
print(result)  # Should confirm Wada property
plot(result, basins = basins)

# Compare boundary complexity at different resolutions
boundary <- get_boundary(basins)
print(nrow(boundary))  # Number of boundary points
} # }
```
