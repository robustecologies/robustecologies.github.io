# Compute normalized solid angle measure in R^3 using Euler-Lagrange formula

Computes the normalized solid angle measure of a cone in R^3 generated
by three unit vectors using the classical Euler-Lagrange formula.

## Usage

``` r
solid_angle_3d(v1, v2, v3)
```

## Arguments

- v1:

  First unit vector (length 3)

- v2:

  Second unit vector (length 3)

- v3:

  Third unit vector (length 3)

## Value

The normalized solid angle measure (between 0 and 1)

## Details

The formula from Euler and Lagrange (page 1 of the paper) is:

\$\$E = 2 \arctan2\left(\|v_1 \cdot (v_2 \times v_3)\|,\\ 1 + v_2 \cdot
v_3 + v_2 \cdot v_1 + v_1 \cdot v_3\right)\$\$

The normalized measure is \\E / (4\pi)\\.

## References

Fitisone, A., & Zhou, Y. (2023). Solid angle measure of polyhedral
cones. arXiv:2304.11102 (math.CO). <https://arxiv.org/abs/2304.11102>

Van Oosterom, A., & Strackee, J. (1983). The solid angle of a plane
triangle. *IEEE Transactions on Biomedical Engineering*, BME-30(2),
125-126.
[doi:10.1109/TBME.1983.325207](https://doi.org/10.1109/TBME.1983.325207)

## Examples

``` r
# Orthogonal vectors (octant = 1/8 of space)
v1 <- c(1, 0, 0)
v2 <- c(0, 1, 0)
v3 <- c(0, 0, 1)
solid_angle_3d(v1, v2, v3)  # Should be 0.125
#> [1] 0.125

# Narrow cone
v1 <- c(1, 0, 0)
v2 <- c(1, 0.1, 0) / sqrt(1.01)
v3 <- c(1, 0, 0.1) / sqrt(1.01)
solid_angle_3d(v1, v2, v3)  # Small value
#> [1] 0.0003959095
```
