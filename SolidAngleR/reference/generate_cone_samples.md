# Vectorised uniform sampling on a spherical cap

Generates multiple uniform random points on the spherical cap of axis
\\\hat{\mu}\\ and half-angle \\\theta_0\\ in one call. Dispatches to the
Rcpp/Armadillo backend `generate_cone_samples_cpp()` for the
inverse-transform method (linear time per sample, low overhead) and
falls back to a pure R loop over
[`generate_cone_sample`](https://robustecologies.github.io/SolidAngleR/reference/generate_cone_sample.md)
for the rejection method or when the C++ backend is unavailable.

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

  Positive integer. Number of samples to generate.

- mu_hat:

  Numeric vector. Central axis direction (will be normalised
  internally).

- theta0:

  Numeric. Maximum planar angle (cone half-angle) in radians, \\0 \<
  \theta_0 \le \pi\\.

- method:

  Character scalar; one of `"inverse"` (default; delegated to the C++
  backend) or `"rejection"` (pure R loop).

## Value

A numeric matrix of dimension `length(mu_hat) %*% n_samples`; each
column is a unit vector on the spherical cap.

## Details

The C++ implementation lives in `src/cone_sampling.cpp` and uses
Armadillo for the linear algebra. It applies the same canonical-axis
construction as the R reference
([`generate_cone_sample`](https://robustecologies.github.io/SolidAngleR/reference/generate_cone_sample.md))
and shares its uniformity guarantees by Arun & Venkatapathi (2025).

## References

Arun, I., & Venkatapathi, M. (2025). An O(n) algorithm for generating
uniform random vectors in n-dimensional cones. *Sankhya A*, 87(2),
327-348.
[doi:10.1007/s13171-025-00387-9](https://doi.org/10.1007/s13171-025-00387-9)

## See also

[`generate_cone_sample`](https://robustecologies.github.io/SolidAngleR/reference/generate_cone_sample.md)
for the single-sample reference;
[`generate_hollow_cone_sample`](https://robustecologies.github.io/SolidAngleR/reference/generate_hollow_cone_sample.md)
for the hollow-cap variant;
[`verify_cone_uniformity`](https://robustecologies.github.io/SolidAngleR/reference/verify_cone_uniformity.md)
for goodness-of-fit on the returned samples;
[`rejection_cost`](https://robustecologies.github.io/SolidAngleR/reference/rejection_cost.md)
for the expected rejection cost in the rejection backend.

## Examples

``` r
samples <- generate_cone_samples(1000, c(0, 0, 1), pi / 4)
dim(samples)
#> [1]    3 1000
```
