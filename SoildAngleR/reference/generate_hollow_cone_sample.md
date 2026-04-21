# Generate uniform random point in hollow cone

Generates a random point uniformly distributed on the surface between
two cones (hollow cone), with central axis \\\hat{\mu}\\ and planar
angles \\\theta_1\\ and \\\theta_2\\.

## Usage

``` r
generate_hollow_cone_sample(
  mu_hat,
  theta1,
  theta2,
  method = c("inverse", "rejection")
)
```

## Arguments

- mu_hat:

  Numeric vector. Central axis direction (will be normalized to unit
  vector).

- theta1:

  Numeric. Inner planar angle in radians, \\0 \le \theta_1 \<
  \theta_2\\.

- theta2:

  Numeric. Outer planar angle in radians, \\0 \< \theta_2 \le \pi\\.

- method:

  Character. Method for generating planar angle: `"inverse"` for inverse
  transform sampling (default), or `"rejection"` for rejection sampling.

## Value

Numeric vector. Unit vector uniformly distributed on the hollow cone
surface.

## Details

The hollow cone is defined as the set of all unit vectors \\\hat{x}\\
satisfying: \$\$\cos\theta_2 \le \hat{x} \cdot \hat{\mu} \le
\cos\theta_1\$\$

The algorithm generates a random angle \\\theta\\ from: \$\$\theta =
\Theta^{-1}(U(\Omega_2 - \Omega_1) + \Omega_1)\$\$ where \\U \sim
\text{Uniform}(0,1)\\, \\\Omega_1 = \Theta(\theta_1)\\, and \\\Omega_2 =
\Theta(\theta_2)\\.

## Examples

``` r
# Generate points in hollow cone (between 30 and 60 degrees)
mu_hat <- c(1, 0, 0)
points <- replicate(1000, 
  generate_hollow_cone_sample(mu_hat, pi/6, pi/3))

# Verify angles are in correct range
angles <- acos(points[1, ])
range(angles)  # Should be approximately [pi/6, pi/3]
#> [1] 0.5261712 1.0466128
hist(angles, main = "Distribution in hollow cone")

```
