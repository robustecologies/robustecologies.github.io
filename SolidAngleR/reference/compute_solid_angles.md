# Vectorised solid-angle computation over a list of cones

Applies
[`compute_solid_angle`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angle.md)
to each element of a list of generator matrices and returns the vector
of normalized solid angles.

## Usage

``` r
compute_solid_angles(cone_list, ...)
```

## Arguments

- cone_list:

  A list whose elements are numeric matrices of cone generators, each
  suitable as the `V` argument of
  [`compute_solid_angle`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angle.md).

- ...:

  Further arguments forwarded to
  [`compute_solid_angle`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angle.md)
  (`method`, `max_terms`, `tol`, `normalize`).

## Value

A numeric vector of normalized solid angles, of length equal to
`length(cone_list)`.

## Details

The wrapper is a thin [`sapply()`](https://rdrr.io/r/base/lapply.html)
over `cone_list`; each element is processed independently with the same
`method` and tolerance settings. The cost is the sum of the per-cone
costs and is dominated by the chosen backend.

## References

Fitisone, A., & Zhou, Y. (2023). Solid angle measure of polyhedral
cones. arXiv:2304.11102 (math.CO). <https://arxiv.org/abs/2304.11102>

## See also

[`compute_solid_angle`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angle.md)
for the single-cone dispatcher;
[`diagnose_cone`](https://robustecologies.github.io/SolidAngleR/reference/diagnose_cone.md)
to inspect spectral conditioning before a large batch run.

## Examples

``` r
# Two 3D cones: the orthant and a non-orthogonal triple
cone1 <- diag(3)
cone2 <- cbind(c(1, 0, 0),
               c(0, 1, 0),
               c(1, 1, 0) / sqrt(2))
compute_solid_angles(list(cone1, cone2))
#> Warning: ⚠ Vectors are not linearly independent - solid angle may be 0
#> [1] 0.125 0.000
```
