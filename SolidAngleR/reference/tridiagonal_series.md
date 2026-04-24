# Solid angle via the tridiagonal Fitisone-Zhou series

Computes the normalized solid angle of a simplicial cone whose Gram
matrix \\V^\top V\\ is tridiagonal, using the simplified series of
Fitisone & Zhou (2023, Theorem 4.1, equation 23). The reduction from
\\\binom{n}{2}\\ to \\n - 1\\ multi-index variables makes this backend
asymptotically faster than the general
[`hypergeometric_series`](https://robustecologies.github.io/SolidAngleR/reference/hypergeometric_series.md)
when the structure is available.

## Usage

``` r
tridiagonal_series(V, max_terms = 1000, tol = 1e-10)
```

## Arguments

- V:

  A square \\n \times n\\ numeric matrix whose columns are the cone
  generators. Need not be unit-normalized; columns are rescaled
  internally.

- max_terms:

  Integer. Hard cap on the number of series terms accumulated before
  truncation. Default `1000`.

- tol:

  Numeric. Per-degree convergence tolerance. Default `1e-10`.

## Value

A list with components

- solid_angle:

  Normalized solid angle in \\\[0, 1\]\\.

- n_terms:

  Total number of terms accumulated.

- converged:

  Logical, `TRUE` when the per-degree contribution fell below `tol`
  before `max_terms` was reached.

- beta:

  Numeric vector of length \\n - 1\\ with the consecutive dot products
  \\\beta_i = v_i \cdot v\_{i+1}\\.

- is_tridiagonal:

  Logical flag, `TRUE` when \\V^\top V\\ passed the tridiagonal test
  inside the function.

## Details

For tridiagonal \\V^\top V\\ the Ribando series collapses to
\$\$\Omega(C) = \frac{\|\det V\|}{(4\pi)^{n/2}} \sum\_{b \in
\mathbb{N}\_0^{n-1}} \frac{(-2)^{\|b\|}}{b!}\\ \Gamma\\\left(\frac{1 +
b_1}{2}\right) \prod\_{i=2}^{n-1} \Gamma\\\left(\frac{1 + b\_{i-1} +
b_i}{2}\right) \Gamma\\\left(\frac{1 + b\_{n-1}}{2}\right)\\
\beta^{b},\$\$ with \\\beta_i = v_i \cdot v\_{i+1}\\. By Theorem 4.1 of
Fitisone & Zhou (2023), tridiagonality of \\V^\top V\\ implies positive
definiteness of the associated matrix \\M_n(C)\\, so the series always
converges absolutely.

If the input does not satisfy the tridiagonal test, the function emits a
warning but still iterates; the result is then a controlled
approximation rather than the true solid angle and should not be
trusted.

## References

Fitisone, A., & Zhou, Y. (2023). Solid angle measure of polyhedral
cones. arXiv:2304.11102 (math.CO). Theorem 4.1 and equation (23).
<https://arxiv.org/abs/2304.11102>

## See also

[`is_tridiagonal`](https://robustecologies.github.io/SolidAngleR/reference/is_tridiagonal.md)
for the tridiagonal predicate;
[`create_tridiagonal_cone`](https://robustecologies.github.io/SolidAngleR/reference/create_tridiagonal_cone.md)
for a constructor of tridiagonal generators from given consecutive
angles;
[`hypergeometric_series`](https://robustecologies.github.io/SolidAngleR/reference/hypergeometric_series.md),
[`hypergeometric_series_nd`](https://robustecologies.github.io/SolidAngleR/reference/hypergeometric_series_nd.md)
for the general series;
[`solid_angle_decomposition`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_decomposition.md)
for the decomposition route;
[`compute_solid_angle`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angle.md)
for the dispatcher.

## Examples

``` r
if (FALSE) { # \dontrun{
# 4D cone with only consecutive dot products non-zero
v1 <- c(1, 0, 0, 0)
v2 <- c(0.8, 0.6, 0, 0)
v3 <- c(0, 0.6, 0.8, 0)
v4 <- c(0, 0, 0.7, 0.7) / sqrt(0.98)
V  <- cbind(v1, v2, v3, v4)
res <- tridiagonal_series(V, max_terms = 500)
res$solid_angle
res$is_tridiagonal     # TRUE
} # }
```
