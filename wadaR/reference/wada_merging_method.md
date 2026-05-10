# Detect Wada basins using the merging method

Implements the merging method from Daza et al. (2018) for testing
whether basins of attraction have the Wada property. The algorithm tests
whether basin boundaries remain identical when basins are merged, using
the Hausdorff distance to compare boundaries. Uses parallel Rcpp/OpenMP
for high performance.

## Usage

``` r
wada_merging_method(basins, x_grid = NULL, y_grid = NULL, verbose = TRUE)
```

## Arguments

- basins:

  Matrix of basin assignments or `wada_basins` object from
  [`compute_basins`](https://robustecologies.github.io/wadaR/reference/compute_basins.md).

- x_grid:

  Optional numeric vector of x coordinates. Extracted automatically from
  `wada_basins` object.

- y_grid:

  Optional numeric vector of y coordinates. Extracted automatically from
  `wada_basins` object.

- verbose:

  Logical. Print progress messages and diagnostics. Default is TRUE.

## Value

A list of class `wada_merging_result` containing:

- is_wada:

  Logical. TRUE if basins satisfy the Wada criterion.

- max_distance:

  Numeric. Maximum Hausdorff distance among all boundary pairs.

- min_distance:

  Numeric. Minimum non-zero Hausdorff distance.

- relative_diff:

  Numeric. \\(d\_{\max} - d\_{\min}) / d\_{\min}\\.

- relative_max:

  Numeric. \\d\_{\max}\\ / phase space diagonal.

- distance_matrix:

  Numeric matrix. Pairwise Hausdorff distances between all slim
  boundaries.

- boundaries:

  List of length \\N_A\\. Each element is a matrix (n x 2) of boundary
  point coordinates.

- n_attractors:

  Integer. Number of attractors.

- boundary_counts:

  Integer vector. Number of boundary points for each merged basin.

## Details

**Mathematical background:**

Consider a dynamical system with \\N_A\\ basins of attraction \\B_1,
B_2, \ldots, B\_{N_A}\\. The Wada property states that: \$\$\partial B_1
= \partial B_2 = \cdots = \partial B\_{N_A}\$\$

A key insight is that if basins have the Wada property, merging any
subset of basins should not change the common boundary. Specifically,
for each basin \\B_i\\, define the merged basin: \$\$M_i = \bigcup\_{j
\neq i} B_j\$\$

The boundary between \\B_i\\ and \\M_i\\ is called the *slim boundary*.
For Wada basins, all slim boundaries are identical: \$\$\partial(B_1,
M_1) = \partial(B_2, M_2) = \cdots = \partial(B\_{N_A}, M\_{N_A})\$\$

**The merging algorithm:**

1.  For each basin \\B_i\\ (\\i = 1, \ldots, N_A\\):

    - Merge all other basins into \\M_i\\

    - Extract the slim boundary \\\partial_i\\ between \\B_i\\ and
      \\M_i\\

2.  Compute the Hausdorff distance between all pairs of slim boundaries

3.  Analyze the distance matrix to determine Wada property

**Hausdorff distance:**

The Hausdorff distance between point sets \\X\\ and \\Y\\ is: \$\$d_H(X,
Y) = \max\\\sup\_{x \in X} d(x, Y), \sup\_{y \in Y} d(y, X)\\\$\$

where \\d(x, Y) = \inf\_{y \in Y} \\x - y\\\\ is the distance from point
\\x\\ to set \\Y\\.

**Wada criterion:**

The basins are classified as Wada if the relative difference
\\(d\_{\max} - d\_{\min}) / d\_{\min} \< 2\\ and the maximum distance
normalized by phase space size satisfies \\d\_{\max} / L \< 0.05\\.
These thresholds (`2.0` and `0.05`) are empirical defaults: Daza et al.
(2018) propose adaptive thresholds tied to grid resolution, but these
constants reproduce the published Wada/non-Wada classification on the
canonical benchmarks (Newton fractal for \\z^3 - 1\\ at \\256 \times
256\\; forced damped pendulum at \\F = 1.66\\ at \\256 \times 256\\; the
same pendulum at \\F = 2.5\\ as a non-Wada control). The package ships
regression tests that verify these classifications; if you change
resolution or system substantially, run those tests on your case before
trusting the boolean.

## References

Daza, A., Wagemakers, A., & Sanjuan, M. A. F. (2018). Ascertaining when
a basin is Wada: The merging method. *Scientific Reports*, 8, 9954.
[doi:10.1038/s41598-018-28119-0](https://doi.org/10.1038/s41598-018-28119-0)

Daza, A., Wagemakers, A., Sanjuan, M. A. F., & Yorke, J. A. (2015).
Testing for basins of Wada. *Scientific Reports*, 5, 16579.
[doi:10.1038/srep16579](https://doi.org/10.1038/srep16579)

Kennedy, J., & Yorke, J. A. (1991). Basins of Wada. *Physica D:
Nonlinear Phenomena*, 51(1-3), 213-225.
[doi:10.1016/0167-2789(91)90234-Z](https://doi.org/10.1016/0167-2789%2891%2990234-Z)

## See also

[`compute_basins`](https://robustecologies.github.io/wadaR/reference/compute_basins.md)
for computing basins of attraction,
[`wada_grid_method`](https://robustecologies.github.io/wadaR/reference/wada_grid_method.md)
for an alternative Wada test,
[`wada_straddle_method`](https://robustecologies.github.io/wadaR/reference/wada_straddle_method.md)
for the saddle-straddle method,
[`hausdorff_distance`](https://robustecologies.github.io/wadaR/reference/hausdorff_distance.md)
for the distance metric,
[`plot.wada_merging_result`](https://robustecologies.github.io/wadaR/reference/plot.wada_merging_result.md)
for visualization.

## Examples

``` r
if (FALSE) { # \dontrun{
# ===================================================================== #
# Example 1: Wada basins in forced damped pendulum                      #
# ===================================================================== #
pendulum <- forced_damped_pendulum(forcing = 1.66)
basins <- compute_basins(pendulum, c(-pi, pi), c(-3, 3), resolution = 300)

result <- wada_merging_method(basins)
print(result)

# Visualize all slim boundaries overlaid
plot(result)

# Compare specific boundary pairs
plot(result, show_all = FALSE, boundary1 = 1, boundary2 = 2)

# ===================================================================== #
# Example 2: Newton fractal (mathematically known Wada)                 #
# ===================================================================== #
newton_basins <- compute_newton_basins(n_roots = 3, resolution = 400)
result_newton <- wada_merging_method(newton_basins)
print(result_newton)  # All boundaries should be nearly identical

# Check distance matrix
print(result_newton$distance_matrix)

# ===================================================================== #
# Example 3: Henon-Heiles escape basins                                 #
# ===================================================================== #
hh <- henon_heiles_system(energy = 0.2)
hh_basins <- compute_basins(hh, c(-0.4, 0.4), c(-0.4, 0.4),
                            resolution = 300, t_max = 200)
result_hh <- wada_merging_method(hh_basins)
print(result_hh)
plot(result_hh)

# ===================================================================== #
# Example 4: Compare full vs partial Wada                               #
# ===================================================================== #
# Full Wada (F = 1.66)
p1 <- forced_damped_pendulum(forcing = 1.66)
b1 <- compute_basins(p1, c(-pi, pi), c(-3, 3), resolution = 200)
r1 <- wada_merging_method(b1, verbose = FALSE)

# Partial Wada (F = 1.71)
p2 <- forced_damped_pendulum(forcing = 1.71)
b2 <- compute_basins(p2, c(-pi, pi), c(-3, 3), resolution = 200)
r2 <- wada_merging_method(b2, verbose = FALSE)

cat(sprintf("F = 1.66: rel_diff = %.4f, Wada: %s\n",
            r1$relative_diff, r1$is_wada))
cat(sprintf("F = 1.71: rel_diff = %.4f, Wada: %s\n",
            r2$relative_diff, r2$is_wada))
} # }
```
