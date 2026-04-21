# Map planar angle to solid angle fraction

Computes the solid angle fraction \\\Omega\\ corresponding to a planar
angle \\\theta\\ in n dimensions.

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

## See also

[`omega_to_theta`](https://robustecologies.github.io/SolidAngleR/reference/omega_to_theta.md)

## Examples

``` r
# Solid angle for 45 degrees in 3D
theta_to_omega(pi/4, 3)
#> [1] 0.1464466

# Half sphere should give 0.5
theta_to_omega(pi/2, 3)
#> [1] 0.5
```
