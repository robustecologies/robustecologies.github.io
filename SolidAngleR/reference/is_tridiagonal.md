# Test whether a matrix has tridiagonal structure

Returns `TRUE` when every entry \\M\[i, j\]\\ with \\\|i - j\| \> 1\\
has absolute value below `tol`. Used by the dispatcher to activate the
tridiagonal-series backend.

## Usage

``` r
is_tridiagonal(M, tol = 1e-12)
```

## Arguments

- M:

  A numeric matrix.

- tol:

  Numeric. Off-band absolute-value threshold below which an entry is
  treated as zero. Default `1e-12`.

## Value

A single logical: `TRUE` when only the main diagonal and the two
adjacent diagonals carry mass, `FALSE` otherwise.

## Details

The tridiagonality of \\V^\top V\\ is the structural precondition for
the simplified series of Fitisone & Zhou (2023, Theorem 4.1), which
reduces the multi-index from \\\binom{n}{2}\\ to \\n - 1\\ variables.
Theorem 4.1 also establishes that if \\V^\top V\\ is tridiagonal, then
the associated matrix \\M_n(C)\\ is automatically positive definite.

## References

Fitisone, A., & Zhou, Y. (2023). Solid angle measure of polyhedral
cones. arXiv:2304.11102 (math.CO). <https://arxiv.org/abs/2304.11102>

## See also

[`compute_dot_product_matrix`](https://robustecologies.github.io/SolidAngleR/reference/compute_dot_product_matrix.md)
for the natural input;
[`tridiagonal_series`](https://robustecologies.github.io/SolidAngleR/reference/tridiagonal_series.md)
for the downstream solver activated by a positive result;
[`create_tridiagonal_cone`](https://robustecologies.github.io/SolidAngleR/reference/create_tridiagonal_cone.md)
for the constructor of tridiagonal cones;
[`compute_solid_angle`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angle.md)
for the dispatcher.

## Examples

``` r
is_tridiagonal(matrix(c(1, 0.5, 0,
                        0.5, 1, 0.3,
                        0, 0.3, 1), nrow = 3))                    # TRUE
#> [1] TRUE
is_tridiagonal(matrix(c(1, 0.5, 0.1,
                        0.5, 1, 0.3,
                        0.1, 0.3, 1), nrow = 3))                  # FALSE
#> [1] FALSE
```
