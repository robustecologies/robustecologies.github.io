# Compute solid angle of a polyhedral cone

Calculates the normalized solid angle subtended by a polyhedral cone
defined by a list of unit vectors using Van Oosterom & Strackee formula
for spherical triangles.

## Usage

``` r
solid_angle_polyhedral(vertices)
```

## Arguments

- vertices:

  Matrix where each row represents a unit vector defining vertices of
  the polyhedral cone. Must have dimensions n x 3 where n \>= 3

## Value

Normalized solid angle (0 to 1, where 1 = full sphere)

## Details

For a spherical triangle defined by three unit vectors \\\mathbf{v}\_1,
\mathbf{v}\_2, \mathbf{v}\_3\\, the solid angle is computed using the
Van Oosterom & Strackee (1983) formula:

\$\$\tan(\Omega/2) = \frac{\|\mathbf{v}\_1 \cdot (\mathbf{v}\_2 \times
\mathbf{v}\_3)\|}{1 + \mathbf{v}\_1 \cdot \mathbf{v}\_2 + \mathbf{v}\_1
\cdot \mathbf{v}\_3 + \mathbf{v}\_2 \cdot \mathbf{v}\_3}\$\$

For polygons with more than three vertices, the polygon is decomposed
into triangular fans from the first vertex, and the solid angles are
summed.

The input vertices are automatically normalized to unit length. The
function handles both convex and non-convex polyhedral cones, though for
non-convex cases the interpretation depends on the triangulation order.

## References

Van Oosterom, A., & Strackee, J. (1983). The solid angle of a plane
triangle. *IEEE Transactions on Biomedical Engineering*, BME-30(2),
125-126.
[doi:10.1109/TBME.1983.325207](https://doi.org/10.1109/TBME.1983.325207)

Mazonka, O. (2012). Solid angle of conical surfaces, polyhedral cones,
and intersecting spherical caps. arXiv:1205.1396v2.
[doi:10.48550/arXiv.1205.1396](https://doi.org/10.48550/arXiv.1205.1396)

## Examples

``` r
# Equilateral triangle on unit sphere
vertices <- matrix(c(
  1, 0, 0,
  0, 1, 0,
  0, 0, 1
), nrow = 3, byrow = TRUE)
solid_angle_polyhedral(vertices)
#> [1] 0.125

# Square on unit sphere (cube face)
vertices <- matrix(c(
  1, 1, 1,
  -1, 1, 1,
  -1, -1, 1,
  1, -1, 1
), nrow = 4, byrow = TRUE)
vertices <- vertices / sqrt(rowSums(vertices^2))
solid_angle_polyhedral(vertices)
#> [1] 0.1666667

if (FALSE) { # \dontrun{
# Regular octahedron face (equilateral triangle)
theta <- acos(1/3)  # Angle for octahedron vertices
vertices <- matrix(c(
  cos(0), sin(0), 0,
  cos(2*pi/3), sin(2*pi/3), 0,
  cos(4*pi/3), sin(4*pi/3), 0
), nrow = 3, byrow = TRUE)
vertices[,3] <- sqrt(1 - rowSums(vertices[,1:2]^2))
omega <- solid_angle_polyhedral(vertices)
cat("Octahedron face solid angle:", omega, "\n")
cat("Should be 1/8 =", 1/8, "\n")
} # }
```
