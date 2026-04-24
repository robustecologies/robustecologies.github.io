# Signed decomposition of a simplicial cone into PD or low-dimensional pieces

Decomposes a simplicial cone into finitely many signed sub-cones using
the line-based decomposition of Fitisone & Zhou (2023, Theorem 3.3,
second decomposition). Each sub-cone is either lower-dimensional (and so
contributes zero to the solid angle) or has a positive definite
associated matrix (and so falls within the convergent regime of
Ribando's hypergeometric series).

## Usage

``` r
decompose_cone(V, method = "line")
```

## Arguments

- V:

  A square \\n \times n\\ numeric matrix whose columns are linearly
  independent cone generators.

- method:

  Character. Currently only `"line"` (the Theorem 3.3 line-based
  decomposition) is implemented.

## Value

A list with components

- cones:

  A list of \\n \times n\\ numeric matrices, one per sub-cone in the
  decomposition.

- signs:

  Integer vector of \\\pm 1\\, the inclusion-exclusion sign of each
  sub-cone.

- n_cones:

  Length of `cones`.

- all_pd:

  Logical, `TRUE` when every full-rank sub-cone has a positive definite
  associated matrix.

## Details

By Theorem 3.3 of Fitisone & Zhou (2023), any simplicial cone admits a
decomposition into at most \\(n-1)!\\ sub-cones (Corollary 3.4), each
falling into one of two terminal classes: (I) affine dimension strictly
less than \\n\\, contributing zero to the solid angle; or (II) positive
definite associated matrix, so that Ribando's series (Theorem 1.5)
converges and the sub-cone admits exact computation.

The reconstruction identity is \$\$\Omega_n(C) = \sum\_{i} s_i \cdot
\Omega_n(C_i),\$\$ where \\s_i \in \\-1, +1\\\\ are the signs returned
in `signs` and \\C_i\\ are the cones returned in `cones`. The wrapper
[`solid_angle_decomposition`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_decomposition.md)
performs this aggregation.

Algorithmically, the decomposition is recursive. The current
implementation chooses the last column \\w_n\\ as the splitting line,
computes the signs \\\delta_i = \mathrm{sign}(w_i \cdot w_n)\\ for \\i
\< n\\, and constructs the sub-cones using equations 14-17 of the source
paper.

## References

Fitisone, A., & Zhou, Y. (2023). Solid angle measure of polyhedral
cones. arXiv:2304.11102 (math.CO). <https://arxiv.org/abs/2304.11102>

## See also

[`solid_angle_decomposition`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_decomposition.md)
for the aggregating wrapper that returns the scalar solid angle;
[`compute_associated_matrix`](https://robustecologies.github.io/SolidAngleR/reference/compute_associated_matrix.md),
[`is_positive_definite`](https://robustecologies.github.io/SolidAngleR/reference/is_positive_definite.md)
for the spectral tests used during decomposition;
[`compute_solid_angle`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angle.md)
for the dispatcher.

## Examples

``` r
if (FALSE) { # \dontrun{
# Decompose a 3D cone with non-positive-definite associated matrix
v1 <- c(1, 0, 0)
v2 <- c(-0.5, 0.8, 0.3)
v3 <- c(0.2, 0.3, 0.9)
V  <- normalize_vectors(cbind(v1, v2, v3))
decomp <- decompose_cone(V)
decomp$n_cones
decomp$all_pd
} # }
```
