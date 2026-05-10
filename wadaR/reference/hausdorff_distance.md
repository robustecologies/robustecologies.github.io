# Compute Hausdorff distance between two point sets

Computes the Hausdorff distance between two finite point sets in
\\\mathbb{R}^2\\, measuring how far the sets are from being identical.

## Usage

``` r
hausdorff_distance(X, Y)
```

## Arguments

- X:

  Numeric matrix of points (n x 2). Each row is a point \\(x, y)\\.

- Y:

  Numeric matrix of points (m x 2). Each row is a point \\(x, y)\\.

## Value

Numeric. The Hausdorff distance between X and Y. Returns `Inf` if either
set is empty.

## Details

The Hausdorff distance is defined as the maximum of the directed
Hausdorff distances in both directions: \$\$d_H(X, Y) =
\max\\\vec{d}\_H(X, Y), \vec{d}\_H(Y, X)\\\$\$

where the directed Hausdorff distance is: \$\$\vec{d}\_H(X, Y) =
\sup\_{x \in X} \inf\_{y \in Y} \\x - y\\\$\$

Intuitively, \\d_H(X, Y)\\ is the minimum radius \\r\\ such that every
point in \\X\\ is within distance \\r\\ of some point in \\Y\\, and vice
versa.

**Properties:** The Hausdorff distance satisfies \\d_H(X, Y) \geq 0\\
(non-negativity), \\d_H(X, Y) = 0\\ iff \\X = Y\\ (identity of
indiscernibles), \\d_H(X, Y) = d_H(Y, X)\\ (symmetry), and \\d_H(X, Z)
\leq d_H(X, Y) + d_H(Y, Z)\\ (triangle inequality).

**Implementation:**

This function uses an efficient C++ implementation with O(n log m + m
log n) complexity. For very large point sets, consider subsampling.

## References

Hausdorff, F. (1914). *Grundzuge der Mengenlehre*. Leipzig: Veit.

Huttenlocher, D. P., Klanderman, G. A., & Rucklidge, W. J. (1993).
Comparing images using the Hausdorff distance. *IEEE Transactions on
Pattern Analysis and Machine Intelligence*, 15(9), 850-863.
[doi:10.1109/34.232073](https://doi.org/10.1109/34.232073)

## See also

[`wada_merging_method`](https://robustecologies.github.io/wadaR/reference/wada_merging_method.md)
which uses this metric for boundary comparison.

## Examples

``` r
# Simple example: two sets of points
X <- matrix(c(0, 0, 1, 0, 0, 1), ncol = 2, byrow = TRUE)
Y <- matrix(c(0.1, 0, 1.1, 0, 0.1, 1), ncol = 2, byrow = TRUE)
d <- hausdorff_distance(X, Y)
print(d)  # Should be ~0.1
#> [1] 0.1

if (FALSE) { # \dontrun{
# Compare Newton fractal boundaries from different basins
newton_basins <- compute_newton_basins(n_roots = 3, resolution = 300)
result <- wada_merging_method(newton_basins)

# Get individual boundaries
b1 <- result$boundaries[[1]]
b2 <- result$boundaries[[2]]

# Compute Hausdorff distance directly
d12 <- hausdorff_distance(b1, b2)
print(d12)

# For Wada basins, d12 should be small
} # }
```
