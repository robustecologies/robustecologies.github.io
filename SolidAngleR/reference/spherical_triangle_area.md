# Area of the spherical triangle on the unit sphere

Computes the unnormalized solid angle (in steradians) of the spherical
triangle on \\S^2\\ whose vertices are three unit vectors. Equivalent to
the spherical excess of the triangle.

## Usage

``` r
spherical_triangle_area(v1, v2, v3)
```

## Arguments

- v1:

  Numeric vector of length three; first triangle vertex.

- v2:

  Numeric vector of length three; second triangle vertex.

- v3:

  Numeric vector of length three; third triangle vertex.

## Value

A single non-negative numeric value: the area of the spherical triangle
in steradians. The full sphere has total area \\4\pi\\, so the
normalized solid angle is `result / (4 * pi)`.

## Details

The function evaluates the same Van Oosterom-Strackee formula used in
[`solid_angle_3d`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_3d.md)
but returns the unnormalized excess \$\$E = 2\\
\mathrm{atan2}\\\left(\|v_1 \cdot (v_2 \times v_3)\|,\\ 1 + v_1 \cdot
v_2 + v_2 \cdot v_3 + v_3 \cdot v_1\right).\$\$ Mathematically, \\E\\
equals the spherical excess \\(\alpha + \beta + \gamma) - \pi\\ of the
spherical triangle, where \\\alpha, \beta, \gamma\\ are the interior
angles.

L'Huilier's theorem provides an equivalent formulation in terms of arc
lengths; the function
[`lhuilier_angle`](https://robustecologies.github.io/SolidAngleR/reference/lhuilier_angle.md)
implements that alternative.

## References

Van Oosterom, A., & Strackee, J. (1983). The solid angle of a plane
triangle. *IEEE Transactions on Biomedical Engineering*, 30(2), 125-126.
[doi:10.1109/TBME.1983.325207](https://doi.org/10.1109/TBME.1983.325207)

Todhunter, I. (1886). *Spherical Trigonometry*, 5th edition. Macmillan
and Co., London. (Chapter on the spherical excess and L'Huilier's
theorem.)

## See also

[`solid_angle_3d`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_3d.md)
for the normalized solid angle of the same triangle;
[`lhuilier_angle`](https://robustecologies.github.io/SolidAngleR/reference/lhuilier_angle.md)
for the L'Huilier formulation;
[`compute_solid_angle`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angle.md)
for the general dispatcher.

## Examples

``` r
# Octant: spherical triangle area is pi/2 steradians
E <- spherical_triangle_area(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))
E                  # ~ 1.5708
#> [1] 1.570796
E / (4 * pi)       # 0.125 normalized
#> [1] 0.125
```
