# General n-dimensional Ribando series for solid angle

Computes the normalized solid angle of a simplicial cone in arbitrary
dimension by enumerating the multi-indices of the Ribando series via a
recursive weak-composition generator. This is the engine that powers
[`hypergeometric_series`](https://robustecologies.github.io/SolidAngleR/reference/hypergeometric_series.md)
for \\n \geq 4\\.

## Usage

``` r
hypergeometric_series_nd(V, max_terms = 1000, tol = 1e-10)
```

## Arguments

- V:

  A square \\n \times n\\ numeric matrix whose columns are the cone
  generators.

- max_terms:

  Integer. Hard cap on the number of series terms accumulated before
  truncation. Default `1000`. The number of terms at total degree \\d\\
  grows as \\\binom{d + N - 1}{N - 1}\\ with \\N = \binom{n}{2}\\, so
  the cap binds quickly for large \\n\\.

- tol:

  Numeric. Per-degree convergence tolerance. Default `1e-10`.

## Value

A list with components

- solid_angle:

  Normalized solid angle in \\\[0, 1\]\\.

- n_terms:

  Total number of terms accumulated.

- converged:

  Logical flag.

- associated_matrix:

  The associated matrix \\M_n(C)\\.

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

## References

Ribando, J. M. (2006). Measuring solid angles beyond dimension three.
*Discrete & Computational Geometry*, 36(3), 479-487.
[doi:10.1007/s00454-006-1253-4](https://doi.org/10.1007/s00454-006-1253-4)

Aomoto, K. (1977). Analytic structure of Schlafli function. *Nagoya
Mathematical Journal*, 68, 1-16.
<https://projecteuclid.org/euclid.nmj/1118786429>

## See also

[`hypergeometric_series`](https://robustecologies.github.io/SolidAngleR/reference/hypergeometric_series.md)
for the dimension-aware dispatcher that calls this function for \\n \geq
4\\;
[`tridiagonal_series`](https://robustecologies.github.io/SolidAngleR/reference/tridiagonal_series.md)
for the simplified \\(n-1)\\-variable series available when \\V^\top V\\
is tridiagonal;
[`solid_angle_decomposition`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_decomposition.md)
for the decomposition route applicable when the associated matrix is not
positive definite;
[`compute_solid_angle`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angle.md)
for the dispatcher.

## Examples

``` r
if (FALSE) { # \dontrun{
# 4D orthant: should yield 1/16 = 0.0625
hypergeometric_series_nd(diag(4))$solid_angle

# Cross-check the closed-form 3D formula against the general enumerator
V3 <- normalize_vectors(matrix(c(1, 0.1, 0.1,
                                 0.1, 1,   0.1,
                                 0.1, 0.1, 1), nrow = 3))
omega_formula <- solid_angle_3d(V3[, 1], V3[, 2], V3[, 3])
omega_series  <- hypergeometric_series_nd(V3)$solid_angle
c(formula = omega_formula, series = omega_series)
} # }
```
