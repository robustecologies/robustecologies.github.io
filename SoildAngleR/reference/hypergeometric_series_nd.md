# General n-dimensional hypergeometric series for solid angle

Implements Ribando's (2006) hypergeometric series formula for arbitrary
dimensions \\n \ge 2\\. This function uses a recursive generator to
iterate over all "weak compositions" of the degree \\d\\ into \\N =
\binom{n}{2}\\ parts, allowing for the calculation of the solid angle
measure for any simplicial cone with a positive-definite associated
matrix.

## Usage

``` r
hypergeometric_series_nd(V, max_terms = 1000, tol = 1e-10)
```

## Arguments

- V:

  An n x n matrix where columns are the unit vectors generating the
  cone.

- max_terms:

  Maximum number of terms to compute (default: 1000). Note that for
  higher dimensions, the number of terms for a given degree grows
  rapidly.

- tol:

  Convergence tolerance (default: 1e-10).

## Value

A list containing:

- `solid_angle`: The normalized solid angle measure (0 to 1).

- `n_terms`: Total number of terms computed.

- `converged`: Logical indicating if the series converged within
  `max_terms`.

- `associated_matrix`: The associated matrix \\M_n(C)\\.

## Details

**Mathematical Formulation:** The function computes the series:
\$\$T\_{\alpha} = \frac{\|\det V\|}{(4\pi)^{n/2}} \sum\_{\mathbf{a} \in
\mathbb{N}\_0^N} \left\[ \frac{(-2)^{\|\mathbf{a}\|}}{\mathbf{a}!}
\prod\_{i=1}^n \Gamma\left(\frac{1 +
\text{deg}\_i(\mathbf{a})}{2}\right) \right\] \alpha^{\mathbf{a}}\$\$

Here \\N = \binom{n}{2}\\ is the number of pairwise dot products;
\\\mathbf{a} = (a\_{12}, a\_{13}, \dots, a\_{n-1,n})\\ is a multi-index
of exponents; \\\|\mathbf{a}\| = \sum\_{j\<k} a\_{jk}\\ is the total
degree of the term; and \\\text{deg}\_i(\mathbf{a}) = \sum\_{k \ne i}
a\_{ik}\\ is the "degree" of vertex \\i\\ in the multigraph defined by
\\\mathbf{a}\\.

**Algorithmic Strategy:** Unlike hardcoded implementations for \\n=3\\,
this function dynamically generates the multi-indices \\\mathbf{a}\\ for
increasing total degrees \\d = 0, 1, 2, \dots\\. It uses a recursive
backtracking algorithm to find all weak compositions of \\d\\ into \\N\\
parts.

**Performance and Limitations:** Regarding convergence, the series
converges absolutely if and only if the associated matrix \\M_n(C)\\ is
positive definite, which typically means the cone vectors must not be
too close to each other (the cone must not be too narrow); if this
condition is violated, the series may produce meaningless results or
fail to converge. In terms of computational complexity, the number of
terms grows combinatorially with \\n\\ and the degree \\d\\, so for \\n
\ge 5\\, calculating high-degree terms becomes computationally
expensive. Convergence speed varies substantially: for cones that are
nearly orthogonal (dot products close to 0), convergence is very fast,
whereas for cones with high correlation between generators (even if
positive definite), convergence can be slow and may require many terms
to achieve high precision; in such cases, Monte Carlo methods might be
more efficient.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example 1: 4D Orthant (Exact value should be 1/16 = 0.0625)
V4 <- diag(4)
res4 <- hypergeometric_series_nd(V4)
print(res4$solid_angle)

# Example 2: Compare 3D result with closed formula
# We use a nearly orthogonal cone to ensure fast convergence for demonstration
V3 <- matrix(c(1, 0.1, 0.1, 0.1, 1, 0.1, 0.1, 0.1, 1), nrow = 3)
V3 <- normalize_vectors(V3)

# Closed formula
omega_formula <- solid_angle_3d(V3[,1], V3[,2], V3[,3])

# General series method
res3 <- hypergeometric_series_nd(V3)

cat(sprintf("Formula: %.6f\nSeries:  %.6f\n",
            omega_formula, res3$solid_angle))
} # }
```
