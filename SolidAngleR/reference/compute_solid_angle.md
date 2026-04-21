# Compute normalized solid angle measure of a polyhedral cone

Main function to compute the normalized solid angle measure of a
polyhedral cone in R^n. Automatically selects the most appropriate
method based on dimension and cone properties.

## Usage

``` r
compute_solid_angle(
  V,
  method = "auto",
  max_terms = 1000,
  tol = 1e-10,
  normalize = TRUE
)
```

## Arguments

- V:

  A matrix where columns are vectors generating the cone. For a
  simplicial cone, V should be n x n with linearly independent columns.

- method:

  Method to use. Options:

  - `"auto"`: Automatically select best method (default)

  - `"formula"`: Use closed formula (R^2 or R^3 only)

  - `"series"`: Use hypergeometric series (requires PD matrix)

  - `"tridiagonal"`: Use tridiagonal series

  - `"decomposition"`: Use decomposition method

- max_terms:

  Maximum terms for series methods (default: 1000)

- tol:

  Convergence tolerance (default: 1e-10)

- normalize:

  Whether to normalize input vectors (default: TRUE)

## Value

The normalized solid angle measure (a number between 0 and 1)

## Details

The function implements multiple methods from Fitisone & Zhou (2023):

**For n = 2:** Uses the angle formula: \\\Omega = \theta/(2\pi)\\

**For n = 3:** Uses the Euler-Lagrange formula: \$\$E = 2
\arctan\left(\frac{\|v_1 \cdot (v_2 \times v_3)\|}{1 + v_2 \cdot v_3 +
v_2 \cdot v_1 + v_1 \cdot v_3}\right)\$\$ \\\Omega = E/(4\pi)\\

**For n \\\geq\\ 4 with positive definite associated matrix:** Uses
Ribando's hypergeometric series (Theorem 1.5)

**For n \\\geq\\ 4 with tridiagonal structure:** Uses simplified series
with n-1 coordinates (Theorem 4.1)

**For general case:** Uses decomposition method (Theorems 3.3 and
Corollary 3.4)

## References

Fitisone, A., & Zhou, Y. (2023). Solid angle measure of polyhedral
cones. arXiv:2304.11102 (math.CO). <https://arxiv.org/abs/2304.11102>

Ribando, J. M. (2006). Measuring solid angles beyond dimension three.
*Discrete & Computational Geometry*, 36(3), 479-487.
[doi:10.1007/s00454-006-1253-4](https://doi.org/10.1007/s00454-006-1253-4)

## Examples

``` r
# ========================================================================== #
# Example 1: 2D cone ####
# ========================================================================== #

v1 <- c(1, 0)
v2 <- c(1, 1) / sqrt(2)
V <- cbind(v1, v2)
omega <- compute_solid_angle(V)
print(omega)  # 0.125 (45 degrees / 360 degrees)
#> [1] 0.125

# ========================================================================== #
# Example 2: 3D orthogonal cone (octant) ####
# ========================================================================== #

V <- diag(3)
omega <- compute_solid_angle(V)
print(omega)  # 0.125 (1/8 of space)
#> [1] 0.125

# ========================================================================== #
# Example 3: 4D orthogonal cone ####
# ========================================================================== #
if (FALSE) { # \dontrun{
V <- diag(4)
omega <- compute_solid_angle(V)
print(omega)  # 0.0625 (1/16 of space)
} # }

# ========================================================================== #
# Example 4: Tridiagonal cone ####
# ========================================================================== #
if (FALSE) { # \dontrun{
angles <- c(pi/3, pi/3, pi/3)
V <- create_tridiagonal_cone(angles)
omega <- compute_solid_angle(V, method = "tridiagonal")
print(omega)
} # }

# ========================================================================== #
# Example 5: General cone requiring decomposition ####
# ========================================================================== #
if (FALSE) { # \dontrun{
set.seed(123)
V <- matrix(rnorm(9), nrow = 3)
omega <- compute_solid_angle(V, method = "auto")
print(omega)
} # }
```
