# Calculate the angle between two 3D vectors

Computes the angle between two vectors in 3D space using the dot product
formula. The function is numerically stable and handles edge cases near
0 and \\\pi\\.

## Usage

``` r
angle_between(v1, v2)
```

## Arguments

- v1:

  A numeric vector of length 3 (unit vector recommended).

- v2:

  A numeric vector of length 3 (unit vector recommended).

## Value

The angle in radians, in the range \\\[0, \pi\]\\.

## Details

The angle \\\theta\\ between two vectors \\\mathbf{v}\_1\\ and
\\\mathbf{v}\_2\\ is computed using the dot product formula:

\$\$\theta = \arccos\left(\frac{\mathbf{v}\_1 \cdot
\mathbf{v}\_2}{\|\mathbf{v}\_1\| \|\mathbf{v}\_2\|}\right)\$\$

For unit vectors, this simplifies to:

\$\$\theta = \arccos(\mathbf{v}\_1 \cdot \mathbf{v}\_2)\$\$

The function clamps the dot product to \\\[-1, 1\]\\ to prevent
numerical errors that could result in `acos` receiving values outside
its domain.

## Examples

``` r
# Orthogonal vectors (90 degrees)
v1 <- c(1, 0, 0)
v2 <- c(0, 1, 0)
angle_between(v1, v2)  # pi/2 radians
#> [1] 1.570796

# Parallel vectors (0 degrees)
v1 <- c(1, 1, 1) / sqrt(3)
v2 <- c(2, 2, 2) / sqrt(12)
angle_between(v1, v2)  # ~0 radians
#> [1] 0

# Opposite vectors (180 degrees)
v1 <- c(1, 0, 0)
v2 <- c(-1, 0, 0)
angle_between(v1, v2)  # pi radians
#> [1] 3.141593
```
