# Spherical excess of a triangle via L'Huilier's theorem

Computes the solid angle (spherical excess, in steradians) of a
spherical triangle on the unit sphere from its three side arc-lengths
using L'Huilier's theorem. The formulation is numerically stable for
both small and near-degenerate triangles where the Girard direct sum
loses precision.

## Usage

``` r
lhuilier_angle(a, b, c)
```

## Arguments

- a:

  Side length A (great circle arc length in radians).

- b:

  Side length B (great circle arc length in radians).

- c:

  Side length C (great circle arc length in radians).

## Value

The solid angle of the spherical triangle in steradians.

## Details

### Historical background

L'Huilier's theorem, named after Swiss mathematician Simon Antoine Jean
L'Huilier (1750-1840), provides a numerically stable formula for
computing the spherical excess \\E\\ (solid angle) of a spherical
triangle.

### Mathematical formula

Given a spherical triangle with side lengths \\a\\, \\b\\, \\c\\
(measured as great circle arcs in radians), the solid angle \\E\\ is:

\$\$\tan\left(\frac{E}{4}\right) = \sqrt{\tan\left(\frac{s}{2}\right)
\tan\left(\frac{s-a}{2}\right) \tan\left(\frac{s-b}{2}\right)
\tan\left(\frac{s-c}{2}\right)}\$\$

where \\s = (a + b + c)/2\\ is the semi-perimeter.

Solving for \\E\\:

\$\$E = 4 \arctan\left(\sqrt{\tan\left(\frac{s}{2}\right)
\tan\left(\frac{s-a}{2}\right) \tan\left(\frac{s-b}{2}\right)
\tan\left(\frac{s-c}{2}\right)}\right)\$\$

### Girard's theorem

The spherical excess \\E\\ is related to the interior angles \\\alpha,
\beta, \gamma\\ of the spherical triangle by Girard's theorem:

\$\$E = \alpha + \beta + \gamma - \pi\$\$

For a triangle on a sphere of radius \\R\\, the area is \\R^2 E\\. On
the unit sphere (\\R = 1\\), the area equals the solid angle \\E\\.

## References

L'Huilier, S. A. J. (1789). *Exposition elementaire des principes des
calculs superieurs*. Berlin.

Todhunter, I. (1886). *Spherical Trigonometry* (5th ed.). Macmillan.

Van Oosterom, A., & Strackee, J. (1983). The solid angle of a plane
triangle. *IEEE Transactions on Biomedical Engineering*, 30(2), 125-126.
[doi:10.1109/TBME.1983.325207](https://doi.org/10.1109/TBME.1983.325207)

## See also

[`spherical_triangle_area`](https://robustecologies.github.io/SolidAngleR/reference/spherical_triangle_area.md)
for the Van Oosterom-Strackee form of the same quantity computed from
vertex vectors;
[`angle_between`](https://robustecologies.github.io/SolidAngleR/reference/angle_between.md)
for the side-length helper that feeds this function;
[`solid_angle_3d_from_rays`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_3d_from_rays.md)
for the wrapper that triangulates a polyhedral cone in \\\mathbb{R}^3\\
into spherical triangles and aggregates their L'Huilier excesses.

## Examples

``` r
# Equilateral spherical triangle with 60-degree sides
lhuilier_angle(pi/3, pi/3, pi/3)        # ~ 0.551 steradians
#> [1] 0.5512856

# Right spherical octant: side arcs = pi/2 each
lhuilier_angle(pi/2, pi/2, pi/2)        # pi/2 (one eighth of the sphere)
#> [1] 1.570796

# Small triangle: spherical excess approaches the Euclidean area
a <- b <- c <- 0.1
s <- (a + b + c) / 2
c(spherical = lhuilier_angle(a, b, c),
  euclidean = sqrt(s * (s - a) * (s - b) * (s - c)))
#>   spherical   euclidean 
#> 0.004335547 0.004330127 
```
