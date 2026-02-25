# Detect Wada basins using the grid method

Implements the grid refinement method from Daza et al. (2015) for
testing whether basins of attraction have the Wada property. The
algorithm tests whether a third color (basin) can always be found
between any two colors at successively finer resolutions. Uses parallel
Rcpp/OpenMP for performance.

## Usage

``` r
wada_grid_method(
  basins,
  system_func = NULL,
  attractors = NULL,
  x_range = NULL,
  y_range = NULL,
  max_refinements = 20,
  max_points = 2^15,
  tolerance = 0.005,
  verbose = TRUE
)
```

## Arguments

- basins:

  Matrix of basin assignments or `wada_basins` object from
  [`compute_basins`](https://robustecologies.github.io/wadaR/reference/compute_basins.md).

- system_func:

  Function defining the dynamical system (optional, for refinement at
  finer scales). Usually not needed when using `wada_basins`.

- attractors:

  List of attractor specifications (required with system_func).

- x_range:

  Numeric vector of length 2. Range of x coordinates. Extracted
  automatically from `wada_basins` object.

- y_range:

  Numeric vector of length 2. Range of y coordinates. Extracted
  automatically from `wada_basins` object.

- max_refinements:

  Integer. Maximum number of bisection refinement steps for each
  boundary segment. Default is 20.

- max_points:

  Integer. Maximum number of points to evaluate per segment. Default is
  \\2^{15} = 32768\\.

- tolerance:

  Numeric. Convergence tolerance. Basin is classified as Wada if
  \\W\_{N_A} \> 1 - \text{tolerance}\\. Default is 0.005.

- verbose:

  Logical. Print progress messages and diagnostics. Default is TRUE.

## Value

A list of class `wada_grid_result` containing:

- is_wada:

  Logical. TRUE if basin satisfies the Wada criterion.

- W_m:

  Numeric vector of length \\N_A\\. Proportion of boundary boxes in
  boundary of exactly \\m\\ basins.

- G_m:

  List of length \\N_A\\. Each element contains indices of boxes
  bordering exactly \\m\\ basins.

- G_counts:

  Integer vector. Count of boxes in each class.

- n_attractors:

  Integer. Number of attractors \\N_A\\.

- convergence_step:

  Integer. Refinement step at which convergence occurred (NA if
  max_refinements reached).

- boundary_classification:

  Integer matrix. Value at `[i,j]` indicates the number of distinct
  basins in the neighborhood of that box.

- n_boundary:

  Integer. Total number of boundary boxes.

- params:

  List of input parameters.

## Details

**Mathematical background:**

Consider a dynamical system with \\N_A \geq 3\\ attractors. Each initial
condition \\\mathbf{x}\_0\\ belongs to exactly one basin of attraction
\\B_k\\, \\k = 1, \ldots, N_A\\. The basins have the *Wada property* if
every point on the boundary of one basin is simultaneously on the
boundary of all basins: \$\$\partial B_1 = \partial B_2 = \cdots =
\partial B\_{N_A}\$\$

**The grid method algorithm:**

For a grid of \\n \times n\\ boxes, define for each box \\B_i\\: \$\$m_i
= \text{number of distinct basins in the neighborhood of } B_i\$\$

The algorithm proceeds as follows:

1.  Classify each grid box by the number of distinct basins in its
    8-connected neighborhood (Moore neighborhood)

2.  Identify boundary boxes where \\1 \< m_i \< N_A\\

3.  For each boundary box, construct a segment to the nearest box of a
    different basin

4.  Apply bisection refinement to find additional basins along the
    segment

5.  Update \\m_i\\ based on discovered basins

**Classification statistics:**

The algorithm computes: \$\$G_m = \|\\i : m_i = m\\\| \quad
\text{(number of boxes bordering } m \text{ basins)}\$\$ \$\$W_m =
\frac{G_m}{\sum\_{j\>1} G_j} \quad \text{(proportion of boundary
boxes)}\$\$

**Wada criterion:**

The basin is classified as Wada if \\W\_{N_A} \approx 1\\, meaning all
boundary boxes border all \\N_A\\ basins. The implementation uses:
\$\$\text{is\\wada} = (W\_{N_A} \> 1 - \text{tolerance})\$\$

**Computational complexity warning:**

The grid method runtime depends on the number of boundary boxes and the
refinement depth. For high-resolution grids (\>400x400) with many
boundary points, computation can be slow. The parameter
`max_refinements` controls the maximum bisection depth; reducing it
speeds up computation at the cost of accuracy. Recommendations:

- Use `max_refinements = 10` for quick screening

- Use `max_refinements = 20` (default) for accurate results

- For very large grids, consider using
  [`wada_merging_method`](https://robustecologies.github.io/wadaR/reference/wada_merging_method.md)
  first for quick screening

## References

Daza, A., Wagemakers, A., Sanjuan, M. A. F., & Yorke, J. A. (2015).
Testing for basins of Wada. *Scientific Reports*, 5, 16579.
[doi:10.1038/srep16579](https://doi.org/10.1038/srep16579)

Kennedy, J., & Yorke, J. A. (1991). Basins of Wada. *Physica D:
Nonlinear Phenomena*, 51(1-3), 213-225.
[doi:10.1016/0167-2789(91)90234-Z](https://doi.org/10.1016/0167-2789%2891%2990234-Z)

Nusse, H. E., & Yorke, J. A. (1996). Wada basin boundaries and basin
cells. *Physica D: Nonlinear Phenomena*, 90(3), 242-261.
[doi:10.1016/0167-2789(95)00249-9](https://doi.org/10.1016/0167-2789%2895%2900249-9)

Aguirre, J., Viana, R. L., & Sanjuan, M. A. F. (2009). Fractal
structures in nonlinear dynamics. *Reviews of Modern Physics*, 81(1),
333-386.
[doi:10.1103/RevModPhys.81.333](https://doi.org/10.1103/RevModPhys.81.333)

## See also

[`compute_basins`](https://robustecologies.github.io/wadaR/reference/compute_basins.md)
for computing basins of attraction,
[`wada_merging_method`](https://robustecologies.github.io/wadaR/reference/wada_merging_method.md)
for an alternative Wada test,
[`wada_straddle_method`](https://robustecologies.github.io/wadaR/reference/wada_straddle_method.md)
for the saddle-straddle method,
[`plot.wada_grid_result`](https://robustecologies.github.io/wadaR/reference/plot.wada_grid_result.md)
for visualization.

## Examples

``` r
if (FALSE) { # \dontrun{
# ===================================================================== #
# Example 1: Wada basins in forced damped pendulum                      #
# ===================================================================== #
pendulum <- forced_damped_pendulum(forcing = 1.66)
basins <- compute_basins(pendulum, c(-pi, pi), c(-3, 3), resolution = 300)

result <- wada_grid_method(basins)
print(result)

# Visualize boundary classification
plot(result, basins = basins)

# Show only Wada boundary points
plot(result, basins = basins, show_wada_only = TRUE)

# ===================================================================== #
# Example 2: Newton fractal (known Wada)                                #
# ===================================================================== #
newton_basins <- compute_newton_basins(n_roots = 3, resolution = 400)
result_newton <- wada_grid_method(newton_basins)
print(result_newton)  # Should show W_3 close to 1
plot(result_newton, basins = newton_basins)

# ===================================================================== #
# Example 3: Henon-Heiles escape basins                                 #
# ===================================================================== #
hh <- henon_heiles_system(energy = 0.2)
hh_basins <- compute_basins(hh, c(-0.4, 0.4), c(-0.4, 0.4),
                            resolution = 300, t_max = 200)
result_hh <- wada_grid_method(hh_basins)
print(result_hh)

# ===================================================================== #
# Example 4: Compare different forcing values                           #
# ===================================================================== #
# Full Wada (F = 1.66)
p1 <- forced_damped_pendulum(forcing = 1.66)
b1 <- compute_basins(p1, c(-pi, pi), c(-3, 3), resolution = 200)
r1 <- wada_grid_method(b1, verbose = FALSE)

# Partial Wada (F = 1.71)
p2 <- forced_damped_pendulum(forcing = 1.71)
b2 <- compute_basins(p2, c(-pi, pi), c(-3, 3), resolution = 200)
r2 <- wada_grid_method(b2, verbose = FALSE)

cat(sprintf("F = 1.66: W_3 = %.4f (Wada: %s)\n",
            r1$W_m[3], r1$is_wada))
cat(sprintf("F = 1.71: W_3 = %.4f (Wada: %s)\n",
            r2$W_m[3], r2$is_wada))
} # }
```
