# Normalized solid angle in R^3 from a generator matrix

Convenience wrapper around
[`solid_angle_3d`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_3d.md)
that accepts the three generators as the columns of a \\3 \times 3\\
matrix instead of three separate vectors.

## Usage

``` r
solid_angle_3d_det(V)
```

## Arguments

- V:

  A \\3 \times 3\\ numeric matrix whose columns are the cone generators.
  Need not be unit-normalized; columns are rescaled internally.

## Value

A single numeric value in \\\[0, 1\]\\, the normalized solid angle of
the cone spanned by the columns of `V`.

## Details

Equivalent to `solid_angle_3d(V[, 1], V[, 2], V[, 3])` after
[`normalize_vectors`](https://robustecologies.github.io/SolidAngleR/reference/normalize_vectors.md)
has been applied to `V`. Use this form when the cone is already
represented as a matrix, for example in batch loops over a list of
cones.

## References

Van Oosterom, A., & Strackee, J. (1983). The solid angle of a plane
triangle. *IEEE Transactions on Biomedical Engineering*, 30(2), 125-126.
[doi:10.1109/TBME.1983.325207](https://doi.org/10.1109/TBME.1983.325207)

## See also

[`solid_angle_3d`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_3d.md)
for the three-vector form;
[`compute_solid_angle`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angle.md)
for the general dispatcher;
[`normalize_vectors`](https://robustecologies.github.io/SolidAngleR/reference/normalize_vectors.md)
for the column-normalization helper.

## Examples

``` r
V <- cbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))
solid_angle_3d_det(V)   # 0.125
#> [1] 0.125
```
