# Calculate solid angle of spherical triangle using L'Huilier's theorem

Computes the solid angle (spherical excess) of a spherical triangle on
the unit sphere using L'Huilier's theorem. This classical formula from
spherical trigonometry relates the sides of a spherical triangle to its
area (solid angle).

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
triangle. *IEEE Transactions on Biomedical Engineering*, BME-30(2),
125-126.
[doi:10.1109/TBME.1983.325207](https://doi.org/10.1109/TBME.1983.325207)

## Examples

``` r
# Equilateral spherical triangle with 60-degree sides
a <- b <- c <- pi/3
E <- lhuilier_angle(a, b, c)
print(E)  # ~0.551 steradians
#> [1] 0.5512856

# Right spherical triangle (90-90-90 degrees on sides)
a <- b <- c <- pi/2
E <- lhuilier_angle(a, b, c)
print(E)  # pi/2 steradians (1/8 of sphere)
#> [1] 1.570796

# Small triangle (linearization check)
# For small triangles, E ≈ area in Euclidean geometry
a <- b <- c <- 0.1
E <- lhuilier_angle(a, b, c)
# Compare with Heron's formula for planar triangle
s_planar <- (a + b + c) / 2
area_planar <- sqrt(s_planar * (s_planar - a) * (s_planar - b) * (s_planar - c))
print(c(E, area_planar))  # Should be very close
#> [1] 0.004335547 0.004330127
```
