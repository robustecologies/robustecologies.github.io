# Planar angle to solid-angle fraction in n dimensions

Maps a planar (cross-sectional) angle \\\theta \in \[0, \pi\]\\ to the
fraction of the \\(n-1)\\-sphere covered by a spherical cap of opening
\\\theta\\, expressed via the regularized incomplete beta function.
Inverse of
[`omega_to_theta`](https://robustecologies.github.io/SolidAngleR/reference/omega_to_theta.md).

## Usage

``` r
theta_to_omega(theta, n)
```

## Arguments

- theta:

  Numeric. Planar angle in radians, \\0 \le \theta \le \pi\\.

- n:

  Integer. Dimension of the space.

## Value

Numeric. Solid angle fraction \\\Omega \in \[0, 1\]\\.

## Details

The function \\\Theta: \mathbb{R} \to \mathbb{R}\\ maps the planar
cross-sectional angle \\\theta\\ to the rotated solid angle fraction
\\\Omega\\ in n dimensions: \$\$\Theta(\theta) = \begin{cases}
\frac{1}{2} I(\sin^2\theta; \frac{n-1}{2}, \frac{1}{2}) & \theta \in
\[0, \frac{\pi}{2}\] \\ 1 - \frac{1}{2} I(\sin^2\theta; \frac{n-1}{2},
\frac{1}{2}) & \theta \in (\frac{\pi}{2}, \pi\] \end{cases}\$\$ where
\\I(x; \alpha, \beta)\\ is the regularized incomplete beta function.

## References

Arun, I., & Venkatapathi, M. (2025). An O(n) algorithm for generating
uniform random vectors in n-dimensional cones. *Sankhya A*, 87(2),
327-348.
[doi:10.1007/s13171-025-00387-9](https://doi.org/10.1007/s13171-025-00387-9)

## See also

[`omega_to_theta`](https://robustecologies.github.io/SolidAngleR/reference/omega_to_theta.md)
for the inverse map;
[`generate_planar_angle_inverse`](https://robustecologies.github.io/SolidAngleR/reference/generate_planar_angle_inverse.md)
for inverse-transform sampling that uses both maps;
[`generate_cone_sample`](https://robustecologies.github.io/SolidAngleR/reference/generate_cone_sample.md)
for the downstream cone sampler.

## Examples

``` r
theta_to_omega(pi / 4, 3)   # 45-degree cap in 3D
#> [1] 0.1464466
theta_to_omega(pi / 2, 3)   # hemisphere = 0.5
#> [1] 0.5
```
