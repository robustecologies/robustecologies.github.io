# Generate random planar angle using inverse transform sampling

Generates a random planar angle \\\theta\\ from the distribution that
produces uniform sampling on a spherical cap.

## Usage

``` r
generate_planar_angle_inverse(theta0, n)
```

## Arguments

- theta0:

  Numeric. Maximum planar angle (cone half-angle) in radians.

- n:

  Integer. Dimension of the space.

## Value

Numeric. Random angle \\\theta \in \[0, \theta_0\]\\.

## Details

The cumulative distribution function of \\\theta\\ is:
\$\$F\_\theta(\theta) = \frac{\Theta(\theta)}{\Theta(\theta_0)}\$\$ This
function generates \\\theta\\ using inverse transform sampling:

1.  Generate \\U \sim \text{Uniform}(0, \Omega_0)\\ where \\\Omega_0 =
    \Theta(\theta_0)\\

2.  Return \\\theta = \Theta^{-1}(U)\\

## Note

This method requires computing \\\Theta\\ and \\\Theta^{-1}\\, which may
be vulnerable to floating-point underflow for large n and small
\\\theta_0\\. For such cases, consider using
[`generate_planar_angle_rejection`](https://robustecologies.github.io/SolidAngleR/reference/generate_planar_angle_rejection.md).

## See also

[`generate_planar_angle_rejection`](https://robustecologies.github.io/SolidAngleR/reference/generate_planar_angle_rejection.md)

## Examples

``` r
# Generate angles for 30-degree cone in 5D
angles <- replicate(1000, generate_planar_angle_inverse(pi/6, 5))
hist(angles, main = "Distribution of planar angles")

```
