# Compute expected rejection rate for high-dimensional rejection sampling

Calculates the expected number of samples needed to generate one
accepted sample when using naive rejection sampling on the full sphere.

## Usage

``` r
rejection_cost(theta0, n)
```

## Arguments

- theta0:

  Numeric. Maximum planar angle of the cone in radians.

- n:

  Integer. Dimension of the space.

## Value

Numeric. Expected number of samples per acceptance.

## Details

For small planar angles \\\theta_0\\, the expected number of samples is
approximately: \$\$\frac{1}{p} \approx \frac{\sqrt{2\pi
e}(n-1)}{\theta_0^{n-1}}\$\$ This demonstrates the exponential growth
with dimension n, making rejection sampling prohibitively expensive in
high dimensions.

## Examples

``` r
# Compare rejection rates across dimensions
dims <- c(2, 5, 10, 20, 50, 100)
rates <- sapply(dims, function(n) rejection_cost(pi/6, n))
plot(dims, log10(rates), type = "b",
     xlab = "Dimension", ylab = "log10(Expected samples)",
     main = "Cost of rejection sampling")

```
