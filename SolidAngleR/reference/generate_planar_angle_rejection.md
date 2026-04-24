# Random planar angle by one-dimensional rejection sampling

Draws a single random planar angle \\\theta \in \[0, \theta_0\]\\ from
the distribution required to produce uniform sampling on a spherical cap
of half-angle \\\theta_0\\ in dimension \\n\\, using one-dimensional
rejection sampling with log-domain comparisons. Numerically stable in
high dimensions where the inverse-transform variant underflows.

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

## References

Arun, I., & Venkatapathi, M. (2025). An O(n) algorithm for generating
uniform random vectors in n-dimensional cones. *Sankhya A*, 87(2),
327-348.
[doi:10.1007/s13171-025-00387-9](https://doi.org/10.1007/s13171-025-00387-9)

## See also

[`generate_planar_angle_inverse`](https://robustecologies.github.io/SolidAngleR/reference/generate_planar_angle_inverse.md)
for the inverse-transform alternative;
[`rejection_cost`](https://robustecologies.github.io/SolidAngleR/reference/rejection_cost.md)
for the expected number of rejections in the rejection-based cone
sampler;
[`generate_cone_sample`](https://robustecologies.github.io/SolidAngleR/reference/generate_cone_sample.md)
for the downstream cone sampler.

## Examples

``` r
# Narrow cap in dimension 100: rejection variant is preferred
angles <- replicate(1000, generate_planar_angle_rejection(pi / 12, 100))
summary(angles)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.2420  0.2580  0.2599  0.2591  0.2610  0.2618 
```
