# Detect Wada basins using the saddle-straddle method

Implements the saddle-straddle method from Wagemakers et al. (2020) for
testing whether basins of attraction have the Wada property. The
algorithm tracks the chaotic saddle on basin boundaries and verifies
that there is only one saddle (implying a single Wada boundary). Uses
parallel Rcpp/OpenMP for high performance.

## Usage

``` r
wada_straddle_method(
  system_func,
  attractors = NULL,
  x_range,
  y_range,
  n_points = 10000,
  straddle_eps = 1e-08,
  dt = 0.01,
  max_iter = 1000,
  params = NULL,
  verbose = TRUE
)
```

## Arguments

- system_func:

  Function(state, t) defining the dynamical system, OR a system object
  returned by
  [`forced_damped_pendulum`](https://robustecologies.github.io/wadaR/reference/forced_damped_pendulum.md)
  or
  [`henon_heiles_system`](https://robustecologies.github.io/wadaR/reference/henon_heiles_system.md).
  When a system object is passed, parameters and attractors are
  automatically extracted.

- attractors:

  List of attractor specifications. Each element should be a list with
  `center` (numeric vector of length 2) and `radius` (numeric). Not
  required if `system_func` is a system object.

- x_range:

  Numeric vector of length 2. Range of x coordinates.

- y_range:

  Numeric vector of length 2. Range of y coordinates.

- n_points:

  Integer. Number of saddle points to compute for each merged basin
  pair. Default is 10000. Higher values give more accurate saddle
  approximation but increase computation time.

- straddle_eps:

  Numeric. Segment size threshold for bisection refinement. Default is
  \\10^{-8}\\. Smaller values give more precise boundary approximation.

- dt:

  Numeric. Time step for trajectory integration. Default is 0.01.

- max_iter:

  Integer. Maximum iterations for trajectory evolution. Default is 1000.

- params:

  List of system parameters (e.g., `damping`, `forcing` for pendulum).
  Not required if `system_func` is a system object.

- verbose:

  Logical. Print progress messages and diagnostics. Default is TRUE.

## Value

A list of class `wada_straddle_result` containing:

- is_wada:

  Logical. TRUE if basins satisfy the Wada criterion.

- saddles:

  List of length \\N_A\\. Each element is a matrix (n x 2) of computed
  saddle point coordinates.

- distance_matrix:

  Numeric matrix. Pairwise Hausdorff distances between all saddle
  approximations.

- diameters:

  Numeric vector. Diameter of each saddle point set.

- max_distance:

  Numeric. Maximum Hausdorff distance among all pairs.

- n_attractors:

  Integer. Number of attractors.

- n_points:

  Integer vector. Number of points computed for each saddle.

## Details

**Mathematical background:**

For dissipative dynamical systems with \\N_A \geq 3\\ coexisting
attractors, the basin boundaries contain a chaotic saddle \\S\\—an
invariant set where trajectories exhibit chaotic behavior before
eventually escaping to an attractor.

The Wada property implies that there is a single connected chaotic
saddle lying on all basin boundaries simultaneously: \$\$S = S_1 = S_2 =
\cdots = S\_{N_A}\$\$

where \\S_i\\ is the chaotic saddle on the boundary between basin
\\B_i\\ and the merged basin \\M_i = \bigcup\_{j \neq i} B_j\\.

**The saddle-straddle algorithm:**

The algorithm approximates the chaotic saddle by:

1.  Find a segment straddling the boundary (endpoints in different
    basins)

2.  Use bisection refinement until segment length \\\< \varepsilon\\

3.  Record the midpoint (approximately on the saddle)

4.  Iterate both endpoints forward in time (segment expands along
    unstable manifold of the saddle)

5.  When segment length exceeds threshold, refine and repeat

6.  Continue until desired number of saddle points collected

The key insight is that points on the chaotic saddle repel trajectories
exponentially fast along the unstable manifold, causing the straddling
segment to stretch. Repeated bisection converges to the saddle.

**Wada criterion:**

Compute the Hausdorff distance between all pairs of computed saddles.
The basins are Wada if: \$\$d_H(S_i, S_j) \ll \text{diam}(S)\$\$

where \\\text{diam}(S)\\ is the diameter of the saddle set. The
implementation uses threshold \\d_H \< 0.01 \cdot \text{diam}\\.

**Computational complexity warning:**

This is the most computationally expensive Wada detection method. For
each of the \\N_A\\ merged basin pairs, the algorithm must:

- Find initial straddling points (random sampling)

- Integrate two trajectories forward while collecting saddle points

- Perform repeated bisection refinements

The total number of ODE integration steps scales approximately as:
\$\$O(N_A \times n\\points \times T/\Delta t)\$\$

For 3 basins with `n_points=10000`, this can take **several minutes to
hours**. Recommendations:

- Use `n_points = 1000-2000` for exploration

- Increase to `n_points = 5000-10000` for accurate saddle approximation

- Use
  [`wada_grid_method`](https://robustecologies.github.io/wadaR/reference/wada_grid_method.md)
  or
  [`wada_merging_method`](https://robustecologies.github.io/wadaR/reference/wada_merging_method.md)
  for faster Wada detection; use this method when dynamical insight is
  needed

## References

Wagemakers, A., Daza, A., & Sanjuan, M. A. F. (2020). The
saddle-straddle method to test for Wada basins. *Communications in
Nonlinear Science and Numerical Simulation*, 84, 105167.
[doi:10.1016/j.cnsns.2020.105167](https://doi.org/10.1016/j.cnsns.2020.105167)

Nusse, H. E., & Yorke, J. A. (1989). A procedure for finding numerical
trajectories on chaotic saddles. *Physica D: Nonlinear Phenomena*,
36(1-2), 137-156.
[doi:10.1016/0167-2789(89)90253-4](https://doi.org/10.1016/0167-2789%2889%2990253-4)

Kennedy, J., & Yorke, J. A. (1991). Basins of Wada. *Physica D:
Nonlinear Phenomena*, 51(1-3), 213-225.
[doi:10.1016/0167-2789(91)90234-Z](https://doi.org/10.1016/0167-2789%2891%2990234-Z)

Daza, A., Wagemakers, A., Sanjuan, M. A. F., & Yorke, J. A. (2015).
Testing for basins of Wada. *Scientific Reports*, 5, 16579.
[doi:10.1038/srep16579](https://doi.org/10.1038/srep16579)

## See also

[`forced_damped_pendulum`](https://robustecologies.github.io/wadaR/reference/forced_damped_pendulum.md)
and
[`henon_heiles_system`](https://robustecologies.github.io/wadaR/reference/henon_heiles_system.md)
for creating system objects,
[`wada_grid_method`](https://robustecologies.github.io/wadaR/reference/wada_grid_method.md)
and
[`wada_merging_method`](https://robustecologies.github.io/wadaR/reference/wada_merging_method.md)
for alternative Wada tests,
[`plot.wada_straddle_result`](https://robustecologies.github.io/wadaR/reference/plot.wada_straddle_result.md)
for visualization.

## Examples

``` r
if (FALSE) { # \dontrun{
# ===================================================================== #
# Example 1: Wada basins in forced damped pendulum                      #
# ===================================================================== #
pendulum <- forced_damped_pendulum(forcing = 1.66)

# Pass system object directly (parameters auto-extracted)
result <- wada_straddle_method(pendulum,
                               x_range = c(-pi, pi),
                               y_range = c(-3, 3),
                               n_points = 5000)
print(result)

# Visualize computed saddles
plot(result)

# Show saddles overlaid on basins
basins <- compute_basins(pendulum, c(-pi, pi), c(-3, 3), resolution = 300)
plot(result, basins = basins)

# ===================================================================== #
# Example 2: Compare with grid and merging methods                      #
# ===================================================================== #
pendulum <- forced_damped_pendulum(forcing = 1.66)
basins <- compute_basins(pendulum, c(-pi, pi), c(-3, 3), resolution = 300)

# Run all three methods
r_grid <- wada_grid_method(basins, verbose = FALSE)
r_merge <- wada_merging_method(basins, verbose = FALSE)
r_straddle <- wada_straddle_method(pendulum,
                                   x_range = c(-pi, pi),
                                   y_range = c(-3, 3),
                                   n_points = 5000,
                                   verbose = FALSE)

cat("Grid method:     Wada =", r_grid$is_wada, "\n")
cat("Merging method:  Wada =", r_merge$is_wada, "\n")
cat("Straddle method: Wada =", r_straddle$is_wada, "\n")

# ===================================================================== #
# Example 3: Partial Wada basins (F = 1.71)                             #
# ===================================================================== #
pendulum_partial <- forced_damped_pendulum(forcing = 1.71)
result_partial <- wada_straddle_method(pendulum_partial,
                                       x_range = c(-pi, pi),
                                       y_range = c(-3, 3),
                                       n_points = 5000)
print(result_partial)  # May show NOT WADA or PARTIALLY WADA

# Visualize individual saddles in faceted plot
plot(result_partial, overlay = FALSE)

# ===================================================================== #
# Example 4: Custom parameters                                          #
# ===================================================================== #
# Higher precision saddle computation
result_hires <- wada_straddle_method(pendulum,
                                     x_range = c(-pi, pi),
                                     y_range = c(-3, 3),
                                     n_points = 20000,
                                     straddle_eps = 1e-10,
                                     max_iter = 2000)
print(result_hires)
} # }
```
