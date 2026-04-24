# Solid-angle fraction to planar angle in n dimensions

Maps a normalized solid-angle fraction \\\Omega \in \[0, 1\]\\ back to
the corresponding planar opening angle \\\theta \in \[0, \pi\]\\.
Inverse of
[`theta_to_omega`](https://robustecologies.github.io/SolidAngleR/reference/theta_to_omega.md),
evaluated through the inverse regularized incomplete beta.

## Usage

``` r
omega_to_theta(omega, n)
```

## Arguments

- omega:

  Numeric. Solid angle fraction, \\0 \le \Omega \le 1\\.

- n:

  Integer. Dimension of the space.

## Value

Numeric. Planar angle \\\theta\\ in radians.

## Details

The inverse function \\\Theta^{-1}\\ is given by:
\$\$\Theta^{-1}(\Omega) = \begin{cases} \arcsin\sqrt{I^{-1}(2\Omega;
\frac{n-1}{2}, \frac{1}{2})} & \Omega \in \[0, \frac{1}{2}\] \\ \pi -
\arcsin\sqrt{I^{-1}(2(1-\Omega); \frac{n-1}{2}, \frac{1}{2})} & \Omega
\in (\frac{1}{2}, 1\] \end{cases}\$\$ where \\I^{-1}(y; \alpha, \beta)\\
is the inverse of the regularized incomplete beta function.

## References

Arun, I., & Venkatapathi, M. (2025). An O(n) algorithm for generating
uniform random vectors in n-dimensional cones. *Sankhya A*, 87(2),
327-348.
[doi:10.1007/s13171-025-00387-9](https://doi.org/10.1007/s13171-025-00387-9)

## See also

[`theta_to_omega`](https://robustecologies.github.io/SolidAngleR/reference/theta_to_omega.md)
for the forward map;
[`generate_planar_angle_inverse`](https://robustecologies.github.io/SolidAngleR/reference/generate_planar_angle_inverse.md)
for inverse-transform sampling that uses both maps.

## Examples

``` r
omega_to_theta(0.5, 3)         # hemisphere -> pi/2
#> [1] 1.570796

theta <- pi / 4
omega_to_theta(theta_to_omega(theta, 3), 3)   # round-trip recovers theta
#> [1] 0.7853982
```
