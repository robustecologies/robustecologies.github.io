# Map solid angle fraction to planar angle

Computes the planar angle \\\theta\\ corresponding to a solid angle
fraction \\\Omega\\ in n dimensions. This is the inverse of
[`theta_to_omega`](https://robustecologies.github.io/SolidAngleR/reference/theta_to_omega.md).

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

## See also

[`theta_to_omega`](https://robustecologies.github.io/SolidAngleR/reference/theta_to_omega.md)

## Examples

``` r
# Get angle for half sphere
omega_to_theta(0.5, 3)  # Should be pi/2
#> [1] 1.570796

# Verify inverse relationship
theta <- pi/4
omega <- theta_to_omega(theta, 3)
omega_to_theta(omega, 3)  # Should equal theta
#> [1] 0.7853982
```
