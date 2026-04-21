# Compute solid angle for multiple cones

Vectorized computation of solid angles for a list of cones.

## Usage

``` r
compute_solid_angles(cone_list, ...)
```

## Arguments

- cone_list:

  List of matrices, each representing a cone

- ...:

  Additional arguments passed to `compute_solid_angle`

## Value

Numeric vector of solid angles

## Examples

``` r
# Multiple 3D cones
cone1 <- diag(3)  # Octant
cone2 <- cbind(c(1,0,0), c(0,1,0), c(1,1,0)/sqrt(2))

cones <- list(cone1, cone2)
omegas <- compute_solid_angles(cones)
#> Warning: ⚠ Vectors are not linearly independent - solid angle may be 0
print(omegas)  # c(0.125, ...)
#> [1] 0.125 0.000
```
