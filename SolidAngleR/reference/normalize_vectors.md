# Column-wise unit normalisation of a matrix

Rescales each column of `V` to unit Euclidean norm. Columns whose
Euclidean norm falls below \\10^{-15}\\ (effectively the zero vector)
trigger an error rather than silently producing `NaN`.

## Usage

``` r
normalize_vectors(V)
```

## Arguments

- V:

  A numeric matrix.

## Value

A numeric matrix of the same dimensions as `V`, with each column
rescaled so that \\\\V\[, j\]\\\_2 = 1\\.

## Details

Used as a preprocessing step throughout the package. The geometry of a
simplicial cone is invariant under positive rescaling of its generators,
so unit-normalising before the dispatch removes a degree of freedom that
would otherwise leak into the diagnostic numbers (e.g. condition number
of the Gram matrix).

## References

Strang, G. (2016). *Introduction to Linear Algebra*, 5th edition.
Wellesley-Cambridge Press. ISBN 978-0980232776.

## See also

[`compute_associated_matrix`](https://robustecologies.github.io/SolidAngleR/reference/compute_associated_matrix.md)
and
[`compute_dot_product_matrix`](https://robustecologies.github.io/SolidAngleR/reference/compute_dot_product_matrix.md)
for downstream consumers;
[`compute_solid_angle`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angle.md)
for the dispatcher that calls this helper as part of its preprocessing.

## Examples

``` r
V <- cbind(c(3, 4), c(1, 1))
V_unit <- normalize_vectors(V)
apply(V_unit, 2, function(v) sqrt(sum(v^2)))   # both equal 1
#> [1] 1 1
```
