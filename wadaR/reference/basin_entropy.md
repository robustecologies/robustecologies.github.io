# Compute basin entropy

Calculates the basin entropy, a Shannon-based measure of final-state
unpredictability arising from the fine-scale structure of basins of
attraction.

## Usage

``` r
basin_entropy(
  basins,
  box_size = 10,
  log_base = 2,
  boundary_only = FALSE,
  x_grid = NULL,
  y_grid = NULL
)
```

## Arguments

- basins:

  Matrix of basin assignments or `wada_basins` object from
  [`compute_basins`](https://robustecologies.github.io/wadaR/reference/compute_basins.md)
  or
  [`compute_newton_basins`](https://robustecologies.github.io/wadaR/reference/compute_newton_basins.md).

- box_size:

  Integer. Size of the boxes (in grid cells) used to compute local
  probabilities. Default is 10. Larger boxes capture larger-scale
  uncertainty.

- log_base:

  Numeric. Base of the logarithm for entropy calculation. Default is 2
  (entropy in bits). Use `exp(1)` for nats.

- boundary_only:

  Logical. If TRUE, only compute entropy for boxes on the basin
  boundary. Default is FALSE.

- x_grid:

  Numeric vector of x coordinates. Extracted automatically from
  `wada_basins` object.

- y_grid:

  Numeric vector of y coordinates. Extracted automatically from
  `wada_basins` object.

## Value

A list of class `basin_entropy_result` containing:

- S_b:

  Numeric. Total basin entropy.

- S_max:

  Numeric. Maximum possible entropy \\\log_b(N_A)\\.

- S_normalized:

  Numeric. Normalized entropy \\S_b / S\_{max}\\.

- S_boundary:

  Numeric. Boundary basin entropy (average entropy of boxes on the
  boundary). NA if `boundary_only = FALSE`.

- entropy_matrix:

  Matrix of local entropy values for each box.

- n_boxes:

  Integer. Total number of boxes.

- n_boundary_boxes:

  Integer. Number of boxes on the boundary.

- box_size:

  Integer. Box size used.

- n_attractors:

  Integer. Number of attractors.

- log_base:

  Numeric. Logarithm base used.

## Details

**Mathematical background:**

The basin entropy, introduced by Daza et al. (2016), measures the
unpredictability of a dynamical system's final state given uncertainty
in the initial condition. For a grid divided into boxes of size
\\\epsilon\\, the entropy of box \\i\\ is: \$\$S_i = -\sum\_{j=1}^{N_A}
p\_{ij} \log_b(p\_{ij})\$\$

where \\p\_{ij}\\ is the fraction of box \\i\\ belonging to basin \\j\\,
\\N_A\\ is the number of attractors, and \\b\\ is the logarithm base.
The total basin entropy is the average over all boxes: \$\$S_b =
\frac{1}{N\_{box}} \sum\_{i=1}^{N\_{box}} S_i\$\$

**Interpretation:**

- \\S_b = 0\\: All boxes are monochromatic (no uncertainty)

- \\S_b = \log_b(N_A)\\: Maximum uncertainty (uniform distribution in
  all boxes)

For Wada basins, \\S_b\\ approaches \\\log_b(N_A)\\ on the boundary,
because boundary boxes contain all \\N_A\\ basins at fine scales.

**Boundary basin entropy:**

The boundary basin entropy \\S\_{bb}\\ considers only boxes intersecting
the basin boundary. This focuses on the unpredictability at the boundary
structure, which is most relevant for Wada basins. For true Wada basins,
\\S\_{bb}\\ should approach \\\log_b(N_A)\\.

## References

Shannon, C. E. (1948). A mathematical theory of communication. *Bell
System Technical Journal*, 27(3), 379-423.
[doi:10.1002/j.1538-7305.1948.tb01338.x](https://doi.org/10.1002/j.1538-7305.1948.tb01338.x)

Daza, A., Wagemakers, A., Georgeot, B., Guery-Odelin, D., & Sanjuan, M.
A. F. (2016). Basin entropy: A new tool to analyze uncertainty in
dynamical systems. *Scientific Reports*, 6, 31416.
[doi:10.1038/srep31416](https://doi.org/10.1038/srep31416)

Daza, A., Wagemakers, A., & Sanjuan, M. A. F. (2022). Classifying basins
of attraction using the basin entropy. *Chaos, Solitons & Fractals*,
159, 112112.
[doi:10.1016/j.chaos.2022.112112](https://doi.org/10.1016/j.chaos.2022.112112)

Daza, A., Shipley, J. O., Dolan, S. R., & Sanjuan, M. A. F. (2018). Wada
structures in a binary black hole system. *Physical Review D*, 98(8),
084050.
[doi:10.1103/PhysRevD.98.084050](https://doi.org/10.1103/PhysRevD.98.084050)

## See also

[`compute_basins`](https://robustecologies.github.io/wadaR/reference/compute_basins.md)
for computing basins of attraction,
[`wada_grid_method`](https://robustecologies.github.io/wadaR/reference/wada_grid_method.md)
for testing the Wada property,
[`wada_merging_method`](https://robustecologies.github.io/wadaR/reference/wada_merging_method.md)
for an alternative Wada test.

## Examples

``` r
if (FALSE) { # \dontrun{
# ===================================================================== #
# Example 1: Basin entropy for forced damped pendulum                   #
# ===================================================================== #
pendulum <- forced_damped_pendulum(forcing = 1.66)
basins <- compute_basins(pendulum, c(-pi, pi), c(-3, 3), resolution = 300)

# Compute basin entropy with default box size
entropy_result <- basin_entropy(basins)
print(entropy_result)

# Visualize local entropy
plot(entropy_result)

# ===================================================================== #
# Example 2: Compare Wada vs non-Wada systems                           #
# ===================================================================== #
# Wada basins (F = 1.66)
p1 <- forced_damped_pendulum(forcing = 1.66)
b1 <- compute_basins(p1, c(-pi, pi), c(-3, 3), resolution = 200)
e1 <- basin_entropy(b1)

# Partially Wada (F = 1.71)
p2 <- forced_damped_pendulum(forcing = 1.71)
b2 <- compute_basins(p2, c(-pi, pi), c(-3, 3), resolution = 200)
e2 <- basin_entropy(b2)

cat(sprintf("F = 1.66: S_b = %.3f (normalized: %.3f)\n",
            e1$S_b, e1$S_normalized))
cat(sprintf("F = 1.71: S_b = %.3f (normalized: %.3f)\n",
            e2$S_b, e2$S_normalized))

# ===================================================================== #
# Example 3: Effect of box size on entropy                              #
# ===================================================================== #
basins <- compute_basins(forced_damped_pendulum(forcing = 1.66),
                         c(-pi, pi), c(-3, 3), resolution = 400)

for (box_size in c(5, 10, 20, 40)) {
    e <- basin_entropy(basins, box_size = box_size)
    cat(sprintf("Box size %d: S_b = %.3f, S_boundary = %.3f\n",
                box_size, e$S_b, e$S_boundary))
}

# ===================================================================== #
# Example 4: Newton fractal entropy                                     #
# ===================================================================== #
newton <- compute_newton_basins(n_roots = 5, resolution = 400)
e_newton <- basin_entropy(newton, box_size = 8)
print(e_newton)

# For Newton fractals with n roots, S_max = log2(n)
cat(sprintf("Expected S_max = log2(5) = %.3f\n", log2(5)))
cat(sprintf("Computed S_max = %.3f\n", e_newton$S_max))
} # }
```
