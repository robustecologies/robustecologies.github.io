# Generate multiple uniform random samples in cone

Efficiently generates multiple random points uniformly distributed on a
spherical cap.

## Usage

``` r
generate_cone_samples(
  n_samples,
  mu_hat,
  theta0,
  method = c("inverse", "rejection")
)
```

## Arguments

- n_samples:

  Integer. Number of samples to generate.

- mu_hat:

  Numeric vector. Central axis direction.

- theta0:

  Numeric. Maximum planar angle (cone half-angle) in radians.

- method:

  Character. Method for generating planar angles.

## Value

Matrix with `n_samples` columns, each containing a unit vector on the
spherical cap. If the C++ backend is unavailable, the function falls
back to the pure R implementation.

## Examples

``` r
# Generate 1000 samples on 3D cone
samples <- generate_cone_samples(1000, c(0, 0, 1), pi/4)
```
