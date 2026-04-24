# Construct a tridiagonal simplicial cone from consecutive angles

Constructs a set of \\n\\ unit vectors in \\\mathbb{R}^n\\ whose Gram
matrix \\V^\top V\\ is tridiagonal with the prescribed consecutive
cosines on the off-diagonal. Supplies the canonical input to
[`tridiagonal_series`](https://robustecologies.github.io/SolidAngleR/reference/tridiagonal_series.md).

## Usage

``` r
create_tridiagonal_cone(angles)
```

## Arguments

- angles:

  A numeric vector of length \\n - 1\\ containing the angles, in
  radians, between consecutive generators. Each entry must lie strictly
  within \\(0, \pi)\\.

## Value

An \\n \times n\\ numeric matrix `V` whose columns are unit vectors
satisfying \\v_i \cdot v\_{i+1} = \cos(\text{angles}\[i\])\\ and \\v_i
\cdot v_j = 0\\ for \\\|i - j\| \> 1\\.

## Details

The construction builds the desired tridiagonal Gram matrix \\G\_{ij} =
\delta\_{ij} + \cos(\text{angles}\[i\]) (\delta\_{i, j-1} + \delta\_{i,
j+1})\\ and recovers `V` as the Cholesky factor: \\V\\ upper triangular
with \\V^\top V = G\\. If the constructed \\G\\ is not positive definite
(which happens when the consecutive angles are too large to admit a
coherent set of generators), the function stops with an informative
error.

## References

Fitisone, A., & Zhou, Y. (2023). Solid angle measure of polyhedral
cones. arXiv:2304.11102 (math.CO). <https://arxiv.org/abs/2304.11102>

## See also

[`tridiagonal_series`](https://robustecologies.github.io/SolidAngleR/reference/tridiagonal_series.md)
for the natural consumer of the returned matrix;
[`is_tridiagonal`](https://robustecologies.github.io/SolidAngleR/reference/is_tridiagonal.md)
for the predicate it satisfies;
[`compute_dot_product_matrix`](https://robustecologies.github.io/SolidAngleR/reference/compute_dot_product_matrix.md)
for inspecting the resulting Gram structure;
[`compute_solid_angle`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angle.md)
for the dispatcher that auto-selects the tridiagonal backend.

## Examples

``` r
if (FALSE) { # \dontrun{
angles <- rep(pi / 3, 3)             # consecutive 60-degree angles in R^4
V <- create_tridiagonal_cone(angles)
is_tridiagonal(t(V) %*% V)            # TRUE
tridiagonal_series(V)$solid_angle
} # }
```
