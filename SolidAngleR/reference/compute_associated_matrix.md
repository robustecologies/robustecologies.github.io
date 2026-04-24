# Associated matrix M_n(C) of a simplicial cone

Constructs the associated matrix \\M_n(C)\\ of Ribando (2006, Definition
1.4) for a simplicial cone with generator matrix `V`. This is the
discriminant on which the convergence of the hypergeometric series
rests: the series converges absolutely if and only if \\M_n(C)\\ is
positive definite (Ribando 2006, Theorem 1.5).

## Usage

``` r
compute_associated_matrix(V)
```

## Arguments

- V:

  A square \\n \times n\\ numeric matrix whose columns are the cone
  generators (need not be unit-normalized; the associated matrix is
  invariant under positive rescaling of each column).

## Value

A symmetric \\n \times n\\ numeric matrix with ones on the diagonal and
\\-\|v_i \cdot v_j\|\\ for the off-diagonal entries.

## Details

The associated matrix is defined entry-wise as \$\$m\_{ii} = 1, \qquad
m\_{ij} = -\|\langle v_i, v_j \rangle\|/(\\v_i\\\\\\v_j\\), \quad i \neq
j.\$\$ Positive definiteness of this matrix is the necessary and
sufficient condition for the absolute convergence of the Ribando series
(Theorem 1.5). The dispatcher
[`compute_solid_angle`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angle.md)
consults this criterion when selecting between the series, tridiagonal,
and decomposition backends.

## References

Ribando, J. M. (2006). Measuring solid angles beyond dimension three.
*Discrete & Computational Geometry*, 36(3), 479-487.
[doi:10.1007/s00454-006-1253-4](https://doi.org/10.1007/s00454-006-1253-4)

Fitisone, A., & Zhou, Y. (2023). Solid angle measure of polyhedral
cones. arXiv:2304.11102 (math.CO). <https://arxiv.org/abs/2304.11102>

## See also

[`is_positive_definite`](https://robustecologies.github.io/SolidAngleR/reference/is_positive_definite.md)
for the convergence test built on `M`;
[`compute_dot_product_matrix`](https://robustecologies.github.io/SolidAngleR/reference/compute_dot_product_matrix.md)
for the unsigned Gram matrix \\V^\top V\\;
[`is_tridiagonal`](https://robustecologies.github.io/SolidAngleR/reference/is_tridiagonal.md)
for the specialised tridiagonal predicate;
[`diagnose_cone`](https://robustecologies.github.io/SolidAngleR/reference/diagnose_cone.md)
for the diagnostic wrapper;
[`compute_solid_angle`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angle.md)
for the dispatcher.

## Examples

``` r
V <- cbind(c(1, 0), c(1, 1) / sqrt(2))
M <- compute_associated_matrix(V)
eigen(M)$values        # all positive: convergence guaranteed
#> [1] 1.7071068 0.2928932
```
