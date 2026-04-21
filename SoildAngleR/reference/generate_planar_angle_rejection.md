# Generate random planar angle using one-dimensional rejection sampling

Generates a random planar angle \\\theta\\ from the distribution that
produces uniform sampling on a spherical cap, using rejection sampling.

## Usage

``` r
generate_planar_angle_rejection(theta0, n)
```

## Arguments

- theta0:

  Numeric. Maximum planar angle (cone half-angle) in radians.

- n:

  Integer. Dimension of the space.

## Value

Numeric. Random angle \\\theta \in \[0, \theta_0\]\\.

## Details

The probability density function of \\\theta\\ is: \$\$f\_\theta(\theta)
= \frac{s\_{n-1}}{s_n \Theta(\theta_0)} \sin^{n-2}(\theta)\$\$ This
function uses one-dimensional rejection sampling with log
transformations to avoid floating-point underflow. The acceptance
condition is: \$\$\log(U) \< (n-2)\[\log(\sin\theta) -
\log(\sin\min(\theta_0, \pi/2))\]\$\$ where \\U \sim
\text{Uniform}(0,1)\\.

The average number of rejections is approximately linear in n,
maintaining the O(n) complexity of the overall algorithm.

## See also

[`generate_planar_angle_inverse`](https://robustecologies.github.io/SolidAngleR/reference/generate_planar_angle_inverse.md)

## Examples

``` r
# Generate angles for narrow cone in high dimensions
# (more stable than inverse transform method)
angles <- replicate(1000, generate_planar_angle_rejection(pi/12, 100))
hist(angles, main = "Distribution of planar angles")

```
