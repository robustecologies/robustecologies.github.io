# Compute solid angle using decomposition method

Combines recursive decomposition with hypergeometric series to compute
solid angle of any simplicial cone, even when the associated Gram matrix
is not positive-definite.

## Usage

``` r
solid_angle_decomposition(V, max_terms = 500, tol = 1e-10)
```

## Arguments

- V:

  n x n matrix of unit vectors defining the simplicial cone

- max_terms:

  Maximum number of terms for hypergeometric series (default: 500)

- tol:

  Convergence tolerance for series (default: 1e-10)

## Value

The normalized solid angle measure \\\Omega/(4\pi) \in \[0,1\]\\ where
\\\Omega\\ is the solid angle in steradians

## Details

This function implements the decomposition method based on Theorem 3.3
and Corollary 3.4 from Ribando (2006). The method handles the general
case where the Gram matrix \\M = V^T V\\ is not positive-definite.

**Algorithm**:

1.  Check if the cone is 2D or 3D with positive-definite Gram matrix. If
    so, use direct hypergeometric formulas

2.  Otherwise, decompose the cone recursively using Theorem 3.3:
    \$\$\text{cone}(V) = \text{cone}(V\_{n-1}) \cup \bigcup\_{i=1}^{n-1}
    (-1)^i \text{cone}(V\_{\[i\]})\$\$ where \\V\_{n-1}\\ is the first
    \\n-1\\ columns and \\V\_{\[i\]}\\ replaces column \\i\\ with the
    last column

3.  Compute solid angles for each decomposed cone (with PD matrices)
    using hypergeometric series:
    [`solid_angle_2d`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_2d.md)
    or
    [`solid_angle_3d`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_3d.md)

4.  Sum the signed contributions according to inclusion-exclusion
    principle

**Mathematical background**:

For a simplicial cone in \\\mathbb{R}^n\\ defined by vectors \\V =
\[v_1, \ldots, v_n\]\\, the solid angle is the measure of the spherical
region traced on the unit sphere. When the Gram matrix is
positive-definite, hypergeometric series provide exact computation. When
not PD, decomposition reduces the problem to PD subcones.

The recursion terminates when all subcones have PD Gram matrices,
guaranteeing convergence. However, for some random configurations, deep
recursion may cause stack overflow (see package limitations in NEWS.md).

**When to use**: This method is appropriate when the Gram matrix \\M =
V^T V\\ is not positive-definite, when other methods fail or are not
applicable, or when theoretical verification of other methods is
required.

**Performance**: Computational complexity is \\O(2^n)\\ in worst case
due to recursive decomposition, but typically much faster for
well-conditioned cones.

## References

Ribando, J. M. (2006). Measuring solid angles beyond dimension three.
*Discrete & Computational Geometry*, 36(3), 479-487.
[doi:10.1007/s00454-006-1253-4](https://doi.org/10.1007/s00454-006-1253-4)

Fitisone, A., & Zhou, Y. (2023). Solid angle measure of polyhedral
cones. arXiv:2304.11102 (math.CO). <https://arxiv.org/abs/2304.11102>

## See also

[`solid_angle_2d`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_2d.md),
[`solid_angle_3d`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_3d.md),
[`compute_associated_matrix`](https://robustecologies.github.io/SolidAngleR/reference/compute_associated_matrix.md),
[`is_positive_definite`](https://robustecologies.github.io/SolidAngleR/reference/is_positive_definite.md)

## Examples

``` r
# Example 1: 2D cone (always works directly)
v1 <- c(1, 0)
v2 <- c(0.5, sqrt(3)/2)
V2 <- cbind(v1, v2)
omega2 <- solid_angle_decomposition(V2)
cat("2D cone solid angle:", omega2, "\n")
#> 2D cone solid angle: 0.1666667 

# Example 2: 3D cone with positive-definite Gram matrix
V3_pd <- matrix(c(
  1, 0, 0,
  0, 1, 0,
  0, 0, 1
), nrow = 3, ncol = 3, byrow = TRUE)
omega3_pd <- solid_angle_decomposition(V3_pd)
cat("3D orthogonal cone (octant):", omega3_pd, "should be 1/8 =", 1/8, "\n")
#> 3D orthogonal cone (octant): 0.125 should be 1/8 = 0.125 

# Example 3: 3D cone requiring decomposition (non-PD Gram matrix)
set.seed(123)
V3_npd <- matrix(rnorm(9), nrow = 3)
V3_npd <- normalize_vectors(V3_npd)

M <- compute_associated_matrix(V3_npd)
cat("Is Gram matrix PD?", is_positive_definite(M), "\n")
#> Is Gram matrix PD? FALSE 

omega3_npd <- solid_angle_decomposition(V3_npd)
#> ○ Decomposed into 2 cones
cat("3D non-PD cone solid angle:", omega3_npd, "\n")
#> 3D non-PD cone solid angle: 0.0562723 

if (FALSE) { # \dontrun{
# Example 4: Compare decomposition with direct method
V <- matrix(c(
  1, 1, 1,
  -1, 1, 1,
  0, 0, 1
), nrow = 3, byrow = TRUE)
V <- normalize_vectors(V)

# Check if we can use direct method
M <- compute_associated_matrix(V)
if (is_positive_definite(M)) {
  omega_direct <- solid_angle_3d(V[, 1], V[, 2], V[, 3])
  omega_decomp <- solid_angle_decomposition(V)
  cat("Direct method:", omega_direct, "\n")
  cat("Decomposition:", omega_decomp, "\n")
  cat("Difference:", abs(omega_direct - omega_decomp), "\n")
} else {
  cat("Gram matrix not PD, using decomposition\n")
  omega_decomp <- solid_angle_decomposition(V)
  cat("Result:", omega_decomp, "\n")
}
} # }
```
