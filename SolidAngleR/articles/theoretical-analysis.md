# A mathematical introduction to solid angle methods for polyhedral cones with the SolidAngleR library

``` r

library(SolidAngleR)
```

## Introduction

The computation of solid angles for polyhedral cones represents a
fundamental problem in computational geometry with applications spanning
data envelopment analysis (DEA), lattice point enumeration, computer
graphics, and multivariate probability. While closed-form solutions
exist for dimensions two and three, the extension to arbitrary
dimensions requires sophisticated mathematical machinery. This analysis
synthesizes four key methodological approaches: hypergeometric series
expansions [\[1\]](#ref1)[\[2\]](#ref2), geometric methods for conical
surfaces [\[3\]](#ref3), and graph-theoretic techniques using spanning
trees [\[4\]](#ref4). We progress systematically from classical results
to cutting-edge computational methods, maintaining rigorous mathematical
depth while providing pedagogical clarity suitable for graduate-level
study.

  

## Fundamental concepts

  

### Planar angles in two dimensions

The planar angle \\\theta\\ at a point is defined as the ratio of arc
length to radius:

\\\theta = \frac{l}{r}\\

where \\l\\ is the arc length on a circle of radius \\r\\. **The full
circle subtends \\2\pi\\ radians**, establishing the normalization for
angular measure in two dimensions. This definition extends naturally to
higher dimensions through the concept of solid angles.

``` r

#### Quarter circle (90 degrees = \u03C0/2 radians)                       ####
v1 <- c(1, 0)
v2 <- c(0, 1)
V_90 <- cbind(v1, v2)

omega_90 <- solid_angle_2d(v1, v2)
planar_table <- data.frame(
  quantity = c("Normalized measure", "Absolute angle (radians)"),
  value = c(omega_90, omega_90 * 2 * pi),
  expected = c(0.25, pi/2)
)
knitr::kable(
  planar_table,
  digits = c(NA, 4, 4),
  col.names = c("quantity", "value", "expected")
)
```

| quantity                 |  value | expected |
|:-------------------------|-------:|---------:|
| Normalized measure       | 0.2500 |   0.2500 |
| Absolute angle (radians) | 1.5708 |   1.5708 |

For colinear but opposite rays, the cone degenerates to a line and the
normalized measure is **0** (not 0.5).

  

### Solid angles in three dimensions

A solid angle \\\Omega\\ quantifies the three-dimensional angular extent
subtended by an object at a point. Analogous to planar angles, it is
defined as:

\\\Omega = \frac{A}{r^2}\\

where \\A\\ is the surface area projected onto a sphere of radius \\r\\
centered at the observation point. The unit of solid angle is the
**steradian (sr)**, a dimensionless quantity with \\\text{sr} =
\text{rad}^2\\.

**Key properties:** A complete sphere subtends \\\Omega = 4\pi\\
steradians; a hemisphere subtends \\\Omega = 2\pi\\ steradians; an
orthant (one-eighth of sphere) subtends \\\Omega = \frac{\pi}{2}\\
steradians.

``` r

#### Orthant (positive octant)                                           ####
V_octant <- diag(3)

omega_octant <- solid_angle_3d(V_octant[, 1], V_octant[, 2], V_octant[, 3])
octant_table <- data.frame(
  quantity = c("Normalized measure", "Absolute solid angle (sr)"),
  value = c(omega_octant, omega_octant * 4 * pi),
  expected = c(0.125, pi/2)
)
knitr::kable(
  octant_table,
  digits = c(NA, 6, 6),
  col.names = c("quantity", "value", "expected")
)
```

| quantity                  |    value | expected |
|:--------------------------|---------:|---------:|
| Normalized measure        | 0.125000 | 0.125000 |
| Absolute solid angle (sr) | 1.570796 | 1.570796 |

The **differential solid angle** in spherical coordinates \\(\theta,
\varphi)\\ is:

\\d\Omega = \sin\theta \\ d\theta \\ d\varphi\\

where \\\theta \in \[0, \pi\]\\ is the polar angle and \\\varphi \in
\[0, 2\pi\]\\ is the azimuthal angle. For an arbitrary oriented surface
\\S\\, the solid angle is computed via the surface integral:

\\\Omega = \iint_S \frac{\hat{r} \cdot \hat{n}}{r^2} \\ dS\\

where \\\hat{r}\\ is the unit vector from the observation point to the
surface element and \\\hat{n}\\ is the unit normal.

  

### Right circular cone formula

For a right circular cone with apex at the origin and half-apex angle
\\\theta\\, the solid angle has the elegant closed form:

\\\boxed{\Omega\_{\text{cone}} = 2\pi(1 - \cos\theta)}\\

**Derivation:** The cone intersects the unit sphere in a circle at
colatitude \\\theta\\. The solid angle equals the area of the spherical
cap:

\\\Omega = \int_0^{2\pi} \int_0^\theta \sin\theta' \\ d\theta' \\
d\varphi = 2\pi\[-\cos\theta'\]\_0^\theta = 2\pi(1 - \cos\theta)\\

``` r

#### Test circular cone formula                                          ####
theta_values <- c(pi/6, pi/4, pi/3, pi/2)
theta_degrees <- theta_values * 180 / pi

cone_table <- data.frame(
  theta_deg = theta_degrees,
  omega_normalized = solid_angle_cone(theta_values),
  formula = (1 - cos(theta_values)) / 2
)
knitr::kable(
  cone_table,
  digits = c(1, 6, 6),
  col.names = c("theta (deg)", "omega (normalized)", "formula")
)
```

| theta (deg) | omega (normalized) |  formula |
|------------:|-------------------:|---------:|
|          30 |           0.066987 | 0.066987 |
|          45 |           0.146447 | 0.146447 |
|          60 |           0.250000 | 0.250000 |
|          90 |           0.500000 | 0.500000 |

For small angles, \\\Omega \approx \pi\theta^2\\, recovering the
intuitive result that solid angle scales as the square of the angular
extent.

  

## Spherical triangles and classical theorems

  

### Spherical triangle geometry

A **spherical triangle** is formed by three great circular arcs on a
sphere’s surface. Unlike planar triangles, spherical triangles exhibit
**positive spherical excess**: the sum of interior angles exceeds
\\\pi\\.

Let \\\alpha, \beta, \gamma\\ denote the interior angles and \\a, b, c\\
the opposite side lengths (measured as arc lengths). The **spherical
excess** \\E\\ is defined as:

\\E = \alpha + \beta + \gamma - \pi\\

  

### Girard’s theorem

**Theorem (Girard, 1626 [\[9\]](#ref9)):** The area \\A\\ of a spherical
triangle on a sphere of radius \\R\\ equals \\R^2\\ times its spherical
excess:

\\\boxed{A = R^2(\alpha + \beta + \gamma - \pi) = R^2 E}\\

For the unit sphere (\\R = 1\\), **area equals excess directly**. This
profound result connects the intrinsic geometry (angles) to extrinsic
geometry (area).

**Proof sketch:** Consider the three lunes (double wedges) formed by
extending each side to a complete great circle. Each lune with angle
\\\psi\\ has area \\4\psi R^2\\. The three lunes cover the sphere once
but cover the triangle and its antipodal image three times each:

\\4\alpha R^2 + 4\beta R^2 + 4\gamma R^2 = 4\pi R^2 + 6A\\

Simplifying yields Girard’s formula.

Since solid angle on the unit sphere equals area, we have:

\\\boxed{\Omega\_{\text{triangle}} = \alpha + \beta + \gamma - \pi}\\

This establishes the fundamental connection between solid angles and
spherical geometry.

``` r

#### Right spherical triangle (all angles = 90\u00B0)                     ####
#### Each angle is \u03C0/2, so E = 3(\u03C0/2) - \u03C0 = \u03C0/2                       ####
v1 <- c(1, 0, 0)
v2 <- c(0, 1, 0)
v3 <- c(0, 0, 1)

E_measured <- spherical_triangle_area(v1, v2, v3)
E_girard <- 3 * (pi/2) - pi  # All angles are 90\u00B0
girard_table <- data.frame(
  quantity = c("Measured area (sr)", "Girard's formula (sr)", "Normalized solid angle"),
  value = c(E_measured, E_girard, E_measured / (4 * pi))
)
knitr::kable(
  girard_table,
  digits = c(NA, 6),
  col.names = c("quantity", "value")
)
```

| quantity               |    value |
|:-----------------------|---------:|
| Measured area (sr)     | 1.570796 |
| Girard’s formula (sr)  | 1.570796 |
| Normalized solid angle | 0.125000 |

  

### L’Huilier’s theorem

**Theorem (L’Huilier, 1789 [\[11\]](#ref11)):** For a spherical triangle
with sides \\a, b, c\\ and semiperimeter \\s = \frac{a+b+c}{2}\\, the
spherical excess satisfies:

\\\boxed{\tan\left(\frac{E}{4}\right) =
\sqrt{\tan\left(\frac{s}{2}\right)\tan\left(\frac{s-a}{2}\right)\tan\left(\frac{s-b}{2}\right)\tan\left(\frac{s-c}{2}\right)}}\\

This is the spherical analogue of Heron’s formula, enabling **area
computation from side lengths alone** without computing angles first.
This computational advantage makes L’Huilier’s theorem particularly
valuable for numerical applications.

``` r

#### Equilateral spherical triangle                                      ####
a <- b <- c <- pi/3  # 60-degree sides

E_lhuilier <- lhuilier_angle(a, b, c)
lhuilier_table <- data.frame(
  quantity = c("Spherical excess (sr)", "Normalized solid angle"),
  value = c(E_lhuilier, E_lhuilier / (4 * pi))
)
knitr::kable(
  lhuilier_table,
  digits = c(NA, 6),
  col.names = c("quantity", "value")
)
```

| quantity               |    value |
|:-----------------------|---------:|
| Spherical excess (sr)  | 0.551286 |
| Normalized solid angle | 0.043870 |

``` r


#### For small triangles, compare with Euclidean Heron's formula         ####
a_small <- b_small <- c_small <- 0.1
E_small <- lhuilier_angle(a_small, b_small, c_small)

#### Heron's formula for planar triangle                                 ####
s_planar <- (a_small + b_small + c_small) / 2
area_planar <- sqrt(s_planar * (s_planar - a_small) * (s_planar - b_small) *
                    (s_planar - c_small))

small_triangle_table <- data.frame(
  quantity = c("Spherical (L'Huilier)", "Planar (Heron)", "Relative difference"),
  value = c(E_small, area_planar, abs(E_small - area_planar) / area_planar)
)
knitr::kable(
  small_triangle_table,
  digits = c(NA, 8),
  col.names = c("quantity", "value")
)
```

| quantity              |      value |
|:----------------------|-----------:|
| Spherical (L’Huilier) | 0.00433555 |
| Planar (Heron)        | 0.00433013 |
| Relative difference   | 0.00125167 |

  

### Spherical trigonometry laws

**Spherical law of cosines (for sides):**

\\\cos c = \cos a \cos b + \sin a \sin b \cos \Gamma\\

where \\\Gamma\\ is the angle opposite side \\c\\.

**Spherical law of cosines (for angles):**

\\\cos \Gamma = -\cos A \cos B + \sin A \sin B \cos c\\

**Spherical law of sines:**

\\\frac{\sin a}{\sin A} = \frac{\sin b}{\sin B} = \frac{\sin c}{\sin
C}\\

**Right spherical triangles** (\\\Gamma = \pi/2\\) satisfy the
**spherical Pythagorean theorem**:

\\\cos c = \cos a \cos b\\

These laws, developed extensively by **Euler** [\[10\]](#ref10) and
**Lagrange** in the 18th century using variational methods, form the
foundation for all spherical angle computations.

  

## The Van Oosterom-Strackee formula for tetrahedra

  

### The breakthrough for 3D solid angles

In 1983, Van Oosterom and Strackee [\[5\]](#ref5) published a landmark
formula for the solid angle subtended by a plane triangle, resolving a
longstanding computational challenge.

**Theorem (Van Oosterom-Strackee, 1983 [\[5\]](#ref5)):** Let
\\\vec{r}\_1, \vec{r}\_2, \vec{r}\_3\\ be the position vectors from an
observation point to the vertices of a triangle. The solid angle
subtended is:

\\\boxed{\tan\left(\frac{\Omega}{2}\right) = \frac{\[\vec{r}\_1 \\
\vec{r}\_2 \\ \vec{r}\_3\]}{r_1 r_2 r_3 + (\vec{r}\_1 \cdot
\vec{r}\_2)r_3 + (\vec{r}\_1 \cdot \vec{r}\_3)r_2 + (\vec{r}\_2 \cdot
\vec{r}\_3)r_1}}\\

where \\\[\vec{r}\_1 \\ \vec{r}\_2 \\ \vec{r}\_3\] = \vec{r}\_1 \cdot
(\vec{r}\_2 \times \vec{r}\_3)\\ is the scalar triple product and \\r_i
= \|\vec{r}\_i\|\\.

  

### Mathematical structure and properties

**Sign conventions:** The scalar triple product can be negative
depending on vertex ordering; for consistent orientation, use
\\\|\det(\[\vec{r}\_1, \vec{r}\_2, \vec{r}\_3\])\|\\ when only magnitude
is needed, but account for sign to determine solid angle direction.

**Degeneracies:** When the denominator is negative while the numerator
is positive, add \\\pi\\ to the arctangent result. This handles cases
where the observation point lies in the plane of the triangle.

**Computational efficiency:** This formula reduces solid angle
computation to 3 dot products, 1 cross product, 1 scalar triple product,
and 1 arctangent evaluation. Van Oosterom and Strackee reported a **3×
speedup** compared to numerical integration methods.

  

### Connection to Euler-Lagrange theory

While published in 1983, this result has 18th-century roots. **Euler**
(1755) used calculus of variations to derive spherical trigonometry
formulas, including implicit forms of this relationship. The modern
formulation provides an explicit, computationally efficient expression
of classical theory.

The formula can be derived from L’Huilier’s theorem applied to the
spherical triangle formed by projecting the planar triangle’s vertices
onto the unit sphere, then using the relationships between Cartesian
coordinates and spherical angles.

  

### Extension to polyhedral cones

For a simplicial cone in \\\mathbb{R}^3\\ spanned by unit vectors
\\\vec{v}\_1, \vec{v}\_2, \vec{v}\_3\\ with apex at the origin, we can
directly apply the Van Oosterom-Strackee formula with \\\vec{r}\_i =
\vec{v}\_i\\.

Let \\a = \vec{v}\_2 \cdot \vec{v}\_3\\, \\b = \vec{v}\_1 \cdot
\vec{v}\_3\\, \\c = \vec{v}\_1 \cdot \vec{v}\_2\\, and \\d = \vec{v}\_1
\cdot (\vec{v}\_2 \times \vec{v}\_3)\\. For unit vectors:

\\\tan\left(\frac{\Omega}{2}\right) = \frac{d}{1 + a + b + c}\\

This simplified form for unit vectors is particularly elegant and
computationally efficient.

``` r

#### Test various 3D cones using Van Oosterom-Strackee/Euler-Lagrange    ####
test_cases <- list(
  list(name = "Orthogonal (octant)",
       V = diag(3),
       expected = 0.125),
  list(name = "60 deg between all pairs",
       V = cbind(c(1,0,0), c(0.5,sqrt(3)/2,0), c(0.5,sqrt(3)/6,sqrt(6)/3)),
       expected = NULL),
  list(name = "Narrow cone",
       V = cbind(c(1,0,0), c(0.99,0.1,0), c(0.99,0,0.1)),
       expected = NULL)
)

vo_table <- data.frame(
  case = character(),
  omega_normalized = numeric(),
  expected = numeric(),
  abs_error = numeric(),
  stringsAsFactors = FALSE
)

for (test in test_cases) {
  V <- normalize_vectors(test$V)
  omega <- solid_angle_3d(V[, 1], V[, 2], V[, 3])
  expected_val <- if (is.null(test$expected)) NA_real_ else test$expected
  error_val <- if (is.null(test$expected)) NA_real_ else abs(omega - test$expected)
  vo_table <- rbind(vo_table, data.frame(
    case = test$name,
    omega_normalized = omega,
    expected = expected_val,
    abs_error = error_val
  ))
}

knitr::kable(
  vo_table,
  digits = c(NA, 6, 6, 2),
  col.names = c("case", "omega (normalized)", "expected", "abs. error")
)
```

| case                     | omega (normalized) | expected | abs. error |
|:-------------------------|-------------------:|---------:|-----------:|
| Orthogonal (octant)      |           0.125000 |    0.125 |          0 |
| 60 deg between all pairs |           0.043870 |       NA |         NA |
| Narrow cone              |           0.000404 |       NA |         NA |

  

## Extension to n-dimensional spaces

  

### The challenge beyond dimension three

**No closed-form solution exists for \\n \geq 4\\.** The solid angle of
an \\n\\-dimensional simplicial cone cannot be expressed using
elementary functions for \\n \> 3\\. This fundamental limitation
necessitates alternative approaches including infinite series expansions
(hypergeometric), numerical integration schemes, recursive decomposition
methods, and special function representations.

  

### n-dimensional solid angle definition

The solid angle in \\\mathbb{R}^n\\ measures the \\(n-1)\\-dimensional
content of the cone’s intersection with the unit \\(n-1)\\-sphere,
normalized by the total sphere surface area.

The **total solid angle** of the unit \\(n-1)\\-sphere is:

\\\boxed{\Omega_n = \frac{2\pi^{n/2}}{\Gamma(n/2)}}\\

**Values for specific dimensions:** For \\n=1\\, \\\Omega_1 = 2\\ (two
points); for \\n=2\\, \\\Omega_2 = 2\pi\\ (circle); for \\n=3\\,
\\\Omega_3 = 4\pi\\ (sphere); for \\n=4\\, \\\Omega_4 = 2\pi^2 \approx
19.74\\ (3-sphere); for \\n=5\\, \\\Omega_5 = \frac{8\pi^2}{3} \approx
26.32\\.

The **differential solid angle** in \\n\\ dimensions is:

\\d\Omega_n = \frac{dS\_{n-1}}{r^{n-1}}\\

where \\dS\_{n-1}\\ is the \\(n-1)\\-dimensional surface area element.

  

### Simplicial cones and polyhedral cones

**Definitions:**

A **cone** \\C \subseteq \mathbb{R}^n\\ satisfies: \\\forall x \in C,
\forall \lambda \> 0: \lambda x \in C\\

A **convex cone** additionally satisfies: \\\forall x, y \in C: x + y
\in C\\

A **polyhedral cone** is the conical hull of finitely many vectors:

\\C = \text{cone}(v_1, \ldots, v_m) = \left\\\sum\_{i=1}^m \lambda_i v_i
: \lambda_i \geq 0\right\\\\

A **simplicial cone** in \\\mathbb{R}^n\\ is generated by exactly \\n\\
linearly independent vectors. Every polyhedral cone can be decomposed
into simplicial cones, making simplicial cones the fundamental building
blocks.

  

### Elementary invariances of the solid angle measure

The normalized solid angle is invariant under three elementary
symmetries acting on the generator matrix \\V\\, all of which can be
deduced directly from the rotational invariance of the spherical
measure. First, **orthogonal invariance**: \\\Omega(QV) = \Omega(V)\\
for every \\Q \in O(n)\\, since \\Q\\ is a rigid motion of the ambient
sphere and the \\(n-1)\\-dimensional surface measure is preserved by
isometries. Second, **antipodal symmetry**: \\\Omega(-V) = \Omega(V)\\,
the special case \\Q = -I\\, which means that the cone \\V \cdot
\mathbb{R}^n\_{\geq 0}\\ and its central reflection \\(-V) \cdot
\mathbb{R}^n\_{\geq 0}\\ subtend identical fractions of \\S^{n-1}\\.
Third, **positive column scaling and permutation**: \\\Omega(VDP) =
\Omega(V)\\ whenever \\D\\ is diagonal with strictly positive entries
and \\P\\ is a permutation matrix, because the cone is determined by the
unordered set of rays \\\\\mathbb{R}\_{\geq 0}\\ v_i\\\\, not by the
lengths or labels of the generators.

The antipodal symmetry admits a one-line proof in the Ribando
representation discussed below. The `"mvn"` branch evaluates \\\Omega(V)
= \Pr(W \geq 0)\\ with \\W \sim \mathcal{N}(0,\\ (V^\top V)^{-1})\\.
Substituting \\V \mapsto -V\\ leaves \\V^\top V\\ unchanged, so the
integrand and the integration region of the Genz-Bretz orthant
probability coincide and the value is preserved exactly. The same
conclusion follows in the geometric formulation: the central inversion
\\x \mapsto -x\\ is an isometry of \\S^{n-1}\\, so it preserves
spherical area.

These invariances are useful in applications. When a cone is naturally
specified up to a sign convention, for instance the feasibility cone of
a generalized Lotka-Volterra system, which is generated by the columns
of \\-\boldsymbol{\alpha}\\, the antipodal symmetry guarantees that
passing \\V\\ or \\-V\\ to
[`compute_solid_angle()`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angle.md)
returns the same number; the choice of sign is then governed by which
presentation makes the geometric identity of the cone explicit at the
call site, not by numerical considerations.

``` r

#### Antipodal symmetry: Omega(V) == Omega(-V) up to QMC noise              ####
set.seed(42)
V <- matrix(rnorm(16), 4, 4)
omega_pos <- compute_solid_angle( V, method = "mvn")
omega_neg <- compute_solid_angle(-V, method = "mvn")

knitr::kable(
  data.frame(
    quantity   = c("Omega(V)", "Omega(-V)", "absolute difference"),
    value      = format(c(omega_pos, omega_neg, abs(omega_pos - omega_neg)),
                        digits = 8)
  ),
  caption = "Numerical confirmation of the antipodal symmetry on a random 4D cone.",
  align   = "lr"
)
```

| quantity            |         value |
|:--------------------|--------------:|
| Omega(V)            | 1.9780468e-02 |
| Omega(-V)           | 1.9721233e-02 |
| absolute difference | 5.9234909e-05 |

Numerical confirmation of the antipodal symmetry on a random 4D cone.
{.table}

  

### Gram matrix and inner product structure

For a simplicial cone spanned by \\n\\ vectors \\\\v_1, \ldots, v_n\\\\,
the **Gram matrix** encodes all geometric information:

\\G = \[g\_{ij}\] \quad \text{where} \quad g\_{ij} = v_i \cdot v_j\\

For **unit vectors**, \\G\\ becomes the matrix of mutual inner products
(cosines of angles):

\\G = \[\alpha\_{ij}\] \quad \text{where} \quad \alpha\_{ij} =
\begin{cases} 1 & i = j \\ \cos\theta\_{ij} & i \neq j \end{cases}\\

The Gram matrix is **symmetric** and **positive semi-definite**
(positive definite for linearly independent vectors). Its determinant
equals:

\\\det(G) = \|\det(V)\|^2\\

where \\V\\ is the matrix with columns \\v_1, \ldots, v_n\\.

``` r

#### Example: 4D orthogonal cone                                         ####
V_4d <- diag(4)
G_4d <- t(V_4d) %*% V_4d

knitr::kable(
  G_4d,
  caption = "Gram matrix for 4D orthogonal cone.",
  col.names = paste0("v", seq_len(ncol(G_4d))),
  align = "c"
)
```

| v1  | v2  | v3  | v4  |
|:---:|:---:|:---:|:---:|
|  1  |  0  |  0  |  0  |
|  0  |  1  |  0  |  0  |
|  0  |  0  |  1  |  0  |
|  0  |  0  |  0  |  1  |

Gram matrix for 4D orthogonal cone. {.table}

``` r


gram_4d_table <- data.frame(
  quantity = "Determinant",
  value = sprintf("%.6f", det(G_4d))
)
knitr::kable(gram_4d_table, caption = "Gram determinant (4D orthogonal cone).")
```

| quantity    | value    |
|:------------|:---------|
| Determinant | 1.000000 |

Gram determinant (4D orthogonal cone). {.table}

``` r


#### Example: non-orthogonal 3D cone                                     ####
v1 <- c(1, 0, 0)
v2 <- c(0.5, sqrt(3)/2, 0)  # 60\u00B0 from v1
v3 <- c(0, 0, 1)
V_3d <- cbind(v1, v2, v3)
G_3d <- t(V_3d) %*% V_3d

knitr::kable(
  round(G_3d, 4),
  caption = "Gram matrix for non-orthogonal 3D cone.",
  col.names = paste0("v", seq_len(ncol(G_3d))),
  align = "c"
)
```

|     | v1  | v2  | v3  |
|:----|:---:|:---:|:---:|
| v1  | 1.0 | 0.5 |  0  |
| v2  | 0.5 | 1.0 |  0  |
| v3  | 0.0 | 0.0 |  1  |

Gram matrix for non-orthogonal 3D cone. {.table}

``` r


gram_3d_table <- data.frame(
  quantity = c("Determinant",
               "v1 . v2 (angle = 60 deg)",
               "v1 . v3 (angle = 90 deg)",
               "v2 . v3 (angle = 90 deg)"),
  value = c(sprintf("%.6f", det(G_3d)),
            sprintf("%.4f", G_3d[1, 2]),
            sprintf("%.4f", G_3d[1, 3]),
            sprintf("%.4f", G_3d[2, 3]))
)
knitr::kable(gram_3d_table, caption = "Determinant and pairwise inner products (non-orthogonal 3D cone).")
```

| quantity                 | value    |
|:-------------------------|:---------|
| Determinant              | 0.750000 |
| v1 . v2 (angle = 60 deg) | 0.5000   |
| v1 . v3 (angle = 90 deg) | 0.0000   |
| v2 . v3 (angle = 90 deg) | 0.0000   |

Determinant and pairwise inner products (non-orthogonal 3D cone).
{.table}

  

## The Ribando-Aomoto hypergeometric series method

  

### Historical context

**Kazuhiko Aomoto** (1977) [\[6\]](#ref6) first derived a hypergeometric
series representation for \\n\\-dimensional solid angles in his work on
the analytic structure of the Schläfli function [\[12\]](#ref12). This
profound result connected solid angle computation to the theory of
hyperlogarithmic functions: iterated integrals with logarithmic poles.

**Jason Ribando** (2006) [\[1\]](#ref1) independently rediscovered this
formula while investigating solid angles beyond dimension three,
providing a more accessible derivation and computational perspective.
The convergence properties and practical applicability were rigorously
established, making the method computationally viable.

  

### The hypergeometric series formula

**Theorem (Ribando-Aomoto):** Let \\\Omega \subseteq \mathbb{R}^n\\ be a
simplicial cone spanned by \\n\\ unit vectors \\\\v_1, \ldots, v_n\\\\,
and let \\V\\ be the matrix with these vectors as columns. Define:

- \\\alpha\_{ij} = v_i \cdot v_j\\ for \\i \neq j\\ (with \\\alpha\_{ii}
  = 1\\)
- \\\vec{\alpha} = (\alpha\_{12}, \alpha\_{13}, \ldots, \alpha\_{1n},
  \alpha\_{23}, \ldots, \alpha\_{n-1,n}) \in \mathbb{R}^{\binom{n}{2}}\\
- For multiindex \\\vec{a} = (a\_{12}, \ldots, a\_{n-1,n}) \in
  \mathbb{N}\_0^{\binom{n}{2}}\\: \\\vec{\alpha}^{\vec{a}} = \prod\_{1
  \leq i \< j \leq n} \alpha\_{ij}^{a\_{ij}}\\

Then the solid angle is given by:

\\\boxed{\Omega = \Omega_n \cdot \frac{\|\det(V)\|}{(4\pi)^{n/2}}
\sum\_{\vec{a} \in \mathbb{N}\_0^{\binom{n}{2}}} c\_{\vec{a}} \cdot
\vec{\alpha}^{\vec{a}}}\\

where the **coefficient** \\c\_{\vec{a}}\\ is:

\\c\_{\vec{a}} = \frac{(-2)^{\sum\_{i\<j} a\_{ij}}}{\prod\_{i\<j}
a\_{ij}!} \prod\_{i=1}^{n} \Gamma\left(\frac{1 + \sum\_{m \neq i}
a\_{im}}{2}\right)\\

Here \\\sum\_{m \neq i} a\_{im}\\ denotes the sum of all \\a\_{ij}\\
where \\j \neq i\\ plus all \\a\_{ji}\\ where \\j \neq i\\ (the total
degree at vertex \\i\\ in the graph of the multiindex).

  

### Mathematical structure

**Hypergeometric character:** This is a multivariable Taylor series with
coefficients involving gamma functions. The ratio of successive terms
has the characteristic form of generalized hypergeometric functions.

**Normalization factor:** The term \\\frac{\Omega_n}{(4\pi)^{n/2}}\\
provides dimensional normalization:

\\\frac{\Omega_n}{(4\pi)^{n/2}} =
\frac{2\pi^{n/2}/\Gamma(n/2)}{(4\pi)^{n/2}} = \frac{1}{2^n \pi^{n/2}
\Gamma(n/2)}\\

**Volume factor:** The determinant \\\|\det(V)\|\\ encodes the “index”
or volume scaling of the cone.

  

### Understanding the summation

The sum runs over all non-negative integer multiindices \\\vec{a} \in
\mathbb{N}\_0^{\binom{n}{2}}\\. For dimension \\n\\:

- Number of variables: \\\binom{n}{2} = \frac{n(n-1)}{2}\\
- Each \\a\_{ij} \geq 0\\
- The series is infinite in \\\binom{n}{2}\\ dimensions

**Example for \\n=4\\:**

- 6 variables: \\\alpha\_{12}, \alpha\_{13}, \alpha\_{14}, \alpha\_{23},
  \alpha\_{24}, \alpha\_{34}\\
- Summation over all \\(a\_{12}, a\_{13}, a\_{14}, a\_{23}, a\_{24},
  a\_{34}) \in \mathbb{N}\_0^6\\
- Computational cost: \\O(N^6)\\ for truncation at degree \\N\\

  

### Gamma function products

For each vertex \\i \in \\1, \ldots, n\\\\, the term
\\\Gamma\left(\frac{1 + \sum\_{m \neq i} a\_{im}}{2}\right)\\ encodes
the spherical geometry. Key properties:

\\\Gamma\left(\frac{1}{2}\right) = \sqrt{\pi}\\

\\\Gamma(k + \tfrac{1}{2}) = \frac{(2k)!}{4^k k!}\sqrt{\pi} \quad
\text{for } k \in \mathbb{N}\_0\\

These relationships connect the series coefficients to combinatorial
structures.

  

### Special case: orthonormal frame

When all spanning vectors are orthonormal (\\\alpha\_{ij} = 0\\ for \\i
\neq j\\), only the term with \\\vec{a} = \vec{0}\\ survives:

\\\Omega = \Omega_n \cdot \frac{1}{(4\pi)^{n/2}} \prod\_{i=1}^{n}
\Gamma\left(\frac{1}{2}\right) = \Omega_n \cdot
\frac{\pi^{n/2}}{(4\pi)^{n/2}} = \frac{\Omega_n}{2^n}\\

This gives \\\Omega = \frac{1}{2^n} \times 4\pi = \frac{\pi}{2}\\ for
\\n=3\\, corresponding to one octant—confirming correctness.

  

### Convergence analysis

**Theorem (Ribando, 2006 [\[1\]](#ref1)):** The hypergeometric series
converges when the matrix

\\M = \[m\_{ij}\] \quad \text{where} \quad m\_{ij} = \begin{cases} 1 & i
= j \\ -\|\alpha\_{ij}\| & i \neq j \end{cases}\\

is **positive definite**.

This is the **positive-definite criterion**. Geometrically, it ensures
the cone is “well-formed” and doesn’t span too large a portion of space.

``` r

#### Orthogonal case: M = I (always PD)                                  ####
V_orth <- diag(4)
M_orth <- compute_associated_matrix(V_orth)

knitr::kable(
  M_orth,
  caption = "Associated matrix for 4D orthant.",
  col.names = paste0("v", seq_len(ncol(M_orth))),
  align = "c"
)
```

| v1  | v2  | v3  | v4  |
|:---:|:---:|:---:|:---:|
|  1  |  0  |  0  |  0  |
|  0  |  1  |  0  |  0  |
|  0  |  0  |  1  |  0  |
|  0  |  0  |  0  |  1  |

Associated matrix for 4D orthant. {.table}

``` r


assoc_orth_table <- data.frame(
  quantity = c("Positive definite",
               "Eigenvalues"),
  value = c(as.character(is_positive_definite(M_orth)),
            paste(round(eigen(M_orth, symmetric = TRUE)$values, 4), collapse = ", "))
)
knitr::kable(assoc_orth_table, caption = "Positive-definiteness diagnostic (4D orthant).")
```

| quantity          | value      |
|:------------------|:-----------|
| Positive definite | TRUE       |
| Eigenvalues       | 1, 1, 1, 1 |

Positive-definiteness diagnostic (4D orthant). {.table}

``` r


#### Non-orthogonal case                                                 ####
v1 <- c(1, 0, 0)
v2 <- c(0.5, sqrt(3)/2, 0)  # 60\u00B0 from v1
v3 <- c(0, 0, 1)
V_nonorth <- cbind(v1, v2, v3)
M_nonorth <- compute_associated_matrix(V_nonorth)

knitr::kable(
  round(M_nonorth, 4),
  caption = "Associated matrix for non-orthogonal 3D cone.",
  col.names = paste0("v", seq_len(ncol(M_nonorth))),
  align = "c"
)
```

|  v1  |  v2  | v3  |
|:----:|:----:|:---:|
| 1.0  | -0.5 |  0  |
| -0.5 | 1.0  |  0  |
| 0.0  | 0.0  |  1  |

Associated matrix for non-orthogonal 3D cone. {.table}

``` r


assoc_nonorth_table <- data.frame(
  quantity = c("Positive definite",
               "Eigenvalues"),
  value = c(as.character(is_positive_definite(M_nonorth)),
            paste(round(eigen(M_nonorth, symmetric = TRUE)$values, 4), collapse = ", "))
)
knitr::kable(assoc_nonorth_table, caption = "Positive-definiteness diagnostic (non-orthogonal 3D cone).")
```

| quantity          | value       |
|:------------------|:------------|
| Positive definite | TRUE        |
| Eigenvalues       | 1.5, 1, 0.5 |

Positive-definiteness diagnostic (non-orthogonal 3D cone). {.table}

**Convergence domain:** The series converges within the natural boundary
for solid angles—precisely when the Gram matrix corresponds to a valid
geometric configuration.

**Practical truncation:** For computational purposes, truncate at total
degree \\N\\:

\\\Omega_N = \Omega_n \cdot \frac{\|\det(V)\|}{(4\pi)^{n/2}}
\sum\_{\substack{\vec{a} \in \mathbb{N}\_0^{\binom{n}{2}} \\
\sum\_{i\<j} a\_{ij} \leq N}} c\_{\vec{a}} \cdot
\vec{\alpha}^{\vec{a}}\\

Convergence rate depends on the spectrum of \\M\\. For nearly orthogonal
vectors (small \\\|\alpha\_{ij}\|\\), convergence is rapid.

``` r

#### Orthogonal cones using hypergeometric series                        ####
orth_results <- data.frame(
  n = integer(),
  computed = numeric(),
  expected = numeric(),
  abs_error = numeric(),
  n_terms = integer()
)

for (n in 2:3) {
  V <- diag(n)
  expected <- 1 / 2^n
  result <- hypergeometric_series(V, max_terms = 100)
  orth_results <- rbind(orth_results, data.frame(
    n = n,
    computed = result$solid_angle,
    expected = expected,
    abs_error = abs(result$solid_angle - expected),
    n_terms = result$n_terms
  ))
}

knitr::kable(
  orth_results,
  digits = c(0, 8, 8, 2, 0),
  col.names = c("n", "computed", "expected", "abs. error", "terms used")
)
```

|   n | computed | expected | abs. error | terms used |
|----:|---------:|---------:|-----------:|-----------:|
|   2 |    0.250 |    0.250 |          0 |         12 |
|   3 |    0.125 |    0.125 |          0 |         20 |

``` r


#### For n=4, use the main interface with auto method selection          ####
V4 <- diag(4)
omega4 <- compute_solid_angle(V4, method = "auto", max_terms = 200)
expected4 <- 1 / 16
auto_table <- data.frame(
  n = 4,
  computed = omega4,
  expected = expected4,
  abs_error = abs(omega4 - expected4)
)
knitr::kable(
  auto_table,
  digits = c(0, 8, 8, 2),
  col.names = c("n", "computed", "expected", "abs. error")
)
```

|   n | computed | expected | abs. error |
|----:|---------:|---------:|-----------:|
|   4 |   0.0625 |   0.0625 |          0 |

  

## Fitisone & Zhou’s decomposition methods

  

### The limitation of direct application

The Ribando-Aomoto formula applies directly only when the
**positive-definite criterion** is satisfied. For arbitrary polyhedral
cones, this condition may fail, causing series divergence. **Fitisone
and Zhou (2023)** [\[2\]](#ref2) resolved this limitation through
decomposition methods.

  

### Main contribution: universal decomposition

**Theorem (Fitisone-Zhou, 2023 [\[2\]](#ref2)):** Every full-dimensional
simplicial cone can be decomposed into a finite family of simplicial
cones satisfying the positive-definite criterion.

This establishes that **any polyhedral cone’s solid angle can be
computed** using the hypergeometric series approach after appropriate
decomposition.

  

### Decomposition method 1: general Brion-Vergne decomposition

The first method leverages **Brion-Vergne decomposition theory**,
originally developed for lattice point counting in polytopes.

**Brion’s theorem:** A polytope can be decomposed into tangent cones at
its vertices. Measure-preserving properties allow reconstruction of
global information from local data.

**Application to cones:**

1.  Start with a simplicial cone \\C\\ that violates the
    positive-definite criterion
2.  Decompose \\C\\ into sub-cones \\\\C_1, \ldots, C_k\\\\ using
    Brion-Vergne methods
3.  Each \\C_i\\ satisfies the positive-definite criterion
4.  Compute \\\Omega(C_i)\\ using Ribando-Aomoto formula
5.  Reconstruct: \\\Omega(C) = \sum\_{i=1}^k \epsilon_i \Omega(C_i)\\
    with signs \\\epsilon_i \in \\-1, +1\\\\ from inclusion-exclusion

The decomposition is **algorithmic** and terminates in finite steps.

  

### Decomposition method 2: tridiagonal structure

The second method introduces a **specialized decomposition** yielding
cones with tridiagonal Gram matrices—a major computational optimization.

**Tridiagonal Gram matrix:** A matrix \\G\\ where \\g\_{ij} = 0\\ for
\\\|i-j\| \> 1\\:

\\G = \begin{bmatrix} 1 & \alpha\_{12} & 0 & \cdots & 0 \\ \alpha\_{12}
& 1 & \alpha\_{23} & \cdots & 0 \\ 0 & \alpha\_{23} & 1 & \cdots & 0 \\
\vdots & \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & 0 & \cdots & 1
\end{bmatrix}\\

**Advantages:** The tridiagonal structure offers reduced parameters
(only \\n-1\\ inner products instead of \\\binom{n}{2} =
\frac{n(n-1)}{2}\\); often faster convergence when consecutive
correlations are moderate (fewer terms needed); and lower complexity
(\\O(N \cdot n)\\ instead of \\O(N^{\binom{n}{2}})\\).

**Mechanism:** The tridiagonal structure corresponds to a special
geometric configuration where each vector interacts only with its
immediate neighbors in the ordering.

  

### Convergence for tridiagonal case

**Theorem (Fitisone-Zhou, 2023 [\[2\]](#ref2)):** For cones with
tridiagonal Gram matrices, the hypergeometric series converges within an
explicitly characterized domain. Convergence rates can be established
analytically.

**Practical implication:** Method 2 is typically **computationally
superior** to Method 1 when applicable. Convergence is still slower when
consecutive dot products are large, so tolerance and `max_terms` remain
important tuning parameters.

``` r

#### Create tridiagonal cone with specified consecutive angles           ####
angles <- c(1.5, 1.5, 1.5)  # radians; moderate correlations for fast convergence
V_trid <- create_tridiagonal_cone(angles)

VTV <- t(V_trid) %*% V_trid
knitr::kable(
  round(VTV, 4),
  col.names = paste0("v", seq_len(ncol(VTV))),
  align = "c"
)
```

|   v1   |   v2   |   v3   |   v4   |
|:------:|:------:|:------:|:------:|
| 1.0000 | 0.0707 | 0.0000 | 0.0000 |
| 0.0707 | 1.0000 | 0.0707 | 0.0000 |
| 0.0000 | 0.0707 | 1.0000 | 0.0707 |
| 0.0000 | 0.0000 | 0.0707 | 1.0000 |

``` r


structure_table <- data.frame(
  is_tridiagonal = is_tridiagonal(VTV)
)
knitr::kable(structure_table)
```

| is_tridiagonal |
|:---------------|
| TRUE           |

``` r


#### Compute solid angle using optimized tridiagonal series              ####
result_trid <- tridiagonal_series(V_trid, max_terms = 500, tol = 1e-8)

#### Compare with general hypergeometric series                          ####
result_general <- hypergeometric_series(V_trid, max_terms = 500, tol = 1e-8)

comparison_table <- data.frame(
  method = c("Tridiagonal series", "General series"),
  omega_normalized = c(result_trid$solid_angle, result_general$solid_angle),
  converged = c(result_trid$converged, result_general$converged),
  n_terms = c(result_trid$n_terms, result_general$n_terms),
  abs_diff_vs_tridiag = c(NA_real_,
                          abs(result_trid$solid_angle - result_general$solid_angle))
)
knitr::kable(
  comparison_table,
  digits = c(NA, 6, NA, 0, 2),
  col.names = c("method", "omega (normalized)", "converged",
                "terms used", "abs. diff vs tridiag")
)
```

|  | method | omega (normalized) | converged | terms used | abs. diff vs tridiag |
|:---|:---|---:|:---|---:|---:|
| i | Tridiagonal series | 0.054536 | TRUE | 286 | NA |
|  | General series | 0.054536 | FALSE | 500 | 0 |

``` r

#### Generate random cone (may not have PD associated matrix)            ####
set.seed(42)
V_random <- matrix(rnorm(9), nrow = 3)
V_random <- normalize_vectors(V_random)

M_random <- compute_associated_matrix(V_random)
is_pd <- is_positive_definite(M_random)

#### Compute using automatic method selection                            ####
omega_auto <- compute_solid_angle(V_random, method = "auto")

#### Try decomposition explicitly                                        ####
omega_decomp <- compute_solid_angle(V_random, method = "decomposition")

decomp_rows <- data.frame(
  quantity = c("Associated matrix is PD",
               if (!is_pd) "Eigenvalues",
               "Solid angle (auto)",
               "Solid angle (decomposition)"),
  value = c(as.character(is_pd),
            if (!is_pd) paste(round(eigen(M_random, symmetric = TRUE)$values, 4),
                              collapse = ", "),
            sprintf("%.6f", omega_auto),
            sprintf("%.6f", omega_decomp))
)
knitr::kable(decomp_rows, caption = "Random 3D cone: diagnostic and solid-angle comparison.")
```

| quantity                    | value                   |
|:----------------------------|:------------------------|
| Associated matrix is PD     | FALSE                   |
| Eigenvalues                 | 1.7647, 1.3365, -0.1012 |
| Solid angle (auto)          | 0.034138                |
| Solid angle (decomposition) | 0.034138                |

Random 3D cone: diagnostic and solid-angle comparison. {.table}

  

### Connection to associated matrix

The **associated matrix** \\M\\ (from convergence criterion) relates
directly to the decomposition strategy:

Eigenvalues of \\M\\ determine convergence rate; the condition number
indicates numerical stability; and positive definiteness is the gateway
criterion for applicability.

When \\M\\ has negative eigenvalues, the cone must be decomposed. The
decomposition algorithms construct sub-cones whose associated matrices
are positive definite by design.

  

## Mazonka’s geometric approach

  

### Philosophy: exact closed forms vs. series

While hypergeometric methods require infinite series, **Mazonka (2012)**
[\[3\]](#ref3) developed a **purely geometric approach** yielding exact
closed-form expressions for certain classes of cones. This provides
complementary insight and computational advantages.

  

### General conical surface formula

For a conical surface defined by a closed curve \\\mathbf{s}(l)\\ on the
unit sphere (parameterized by arc length \\l\\), with apex at the
origin:

\\\boxed{\Omega = 2\pi - \sum_i \delta_i - \oint\_{\text{smooth}} dl
\sqrt{\mathbf{u}^2 - (\mathbf{u} \cdot \mathbf{s})^2}}\\

where:

- \\\mathbf{u} = \frac{d\mathbf{s}}{dl}\\ is the unit tangent vector
- \\\delta_i\\ are turn angles at corner vertices (supplementary to
  interior angles)
- The integral is over smooth portions of the curve

**Geometric interpretation:** The formula generalizes the spherical
polygon formula (\\\Omega = \sum \delta_i - (n-2)\pi\\) by adding a line
integral accounting for curvature along smooth arcs.

  

### Polyhedral cones: complex number method

For a polyhedral cone defined by \\n\\ unit vectors \\\\\mathbf{s}\_1,
\ldots, \mathbf{s}\_n\\\\ (in cyclic order), Mazonka derives a
remarkably efficient formula using complex arithmetic.

**Define for each vertex \\j\\:**

\\a_j = \mathbf{s}\_j \cdot \mathbf{s}\_{j-1}, \quad b_j = \mathbf{s}\_j
\cdot \mathbf{s}\_{j+1}\\

\\c_j = \mathbf{s}\_j \cdot (\mathbf{s}\_{j-1} \times
\mathbf{s}\_{j+1})\\

\\d_j = \[\mathbf{s}\_{j-1}, \mathbf{s}\_j, \mathbf{s}\_{j+1}\] =
\mathbf{s}\_{j-1} \cdot (\mathbf{s}\_j \times \mathbf{s}\_{j+1})\\

**Main formula:**

\\\boxed{\Omega = 2\pi - \arg\left\\\prod\_{j=1}^{n} (cb_j - a_j + d_j
i)\right\\}\\

where the product is of complex numbers \\z_j = (cb_j - a_j) + d_j i\\.

**Key advantage:** This requires only **one arctangent operation** (for
the final argument), compared to one per vertex in traditional methods.
Complexity is **O(n)** with constants dominated by vector operations,
not transcendental functions.

  

### Special case: spherical triangles

For \\n=3\\ (triangular face), all \\d_j\\ are equal: \\d_1 = d_2 = d_3
= d\\. This allows factorization:

\\\Omega = 2\pi - \arg\[(bc-a+di)(ac-b+di)(ab-c+di)\]\\

Using \\d^2 = 1 + 2abc - a^2 - b^2 - c^2\\ and simplifying:

\\\Omega = 2\pi - 2\arg\[(1+a+b+c) + di\]\\

Therefore:

\\\boxed{\tan\frac{\Omega}{2} = \frac{d}{1+a+b+c}}\\

**This exactly recovers the Van Oosterom-Strackee formula**,
demonstrating that Mazonka’s approach unifies and extends classical
results.

  

### Intersecting spherical caps

Mazonka addresses a problem not covered by series methods: computing the
solid angle of the **intersection of two cones** with apertures
\\\theta_1, \theta_2\\ and axes separated by angle \\\alpha\\.

**Geometric construction:** Each cone intersects the unit sphere in a
circle (spherical cap boundary); the two circles intersect at two
points; the plane through these intersection points cuts each cap into
segments; the total intersection solid angle equals the sum of the two
cap segments.

**Determining \\\gamma_1\\:** Using spherical trigonometry on the right
triangles formed, the robust form is:

\\\boxed{\gamma_1 = \operatorname{atan2}\\\left(\cos\alpha
\cos\theta_1 - \cos\theta_2,\\ \sin\alpha \cos\theta_1\right)}\\

and analogously for \\\gamma_2\\ (swap indices). The
\\\operatorname{atan2}\\ form preserves sign and quadrant information
that can be lost in a plain \\\tan^{-1}\\ expression.

**Total intersection:**

\\\boxed{\Omega\_{\text{int}} = \Omega(\theta_1, \gamma_1) +
\Omega(\theta_2, \gamma_2)}\\

This provides a closed-form solution for acute cones; practical
implementations also handle special cases (hemispheres, containment,
opposite axes) and may fall back to a robust spherical-cap overlap
integral when the analytic formula is ill-conditioned.

  

## Graph-theoretic approach: spanning trees for extreme rays

  

### Intersection form vs. sum form

**Malekmohammadi and Mostafaee (2017)** [\[4\]](#ref4) introduce a
fundamentally different perspective: converting between cone
representations.

**Intersection form (H-representation):**

\\C\_{\text{int}} = \\x \in \mathbb{R}^n : A_i^T x \leq 0, \\
i=1,\ldots,m\\\\

The cone is defined as the intersection of \\m\\ half-spaces.

**Sum form (V-representation):**

\\C\_{\text{sum}} = \left\\\sum\_{j=1}^k \lambda_j r_j : \lambda_j \geq
0\right\\\\

The cone is the conical hull of its **extreme rays** \\\\r_1, \ldots,
r_k\\\\.

**Fundamental problem:** Given \\C\_{\text{int}}\\, find all extreme
rays \\\\r_j\\\\ to obtain \\C\_{\text{sum}}\\.

  

### Extreme ray characterization

**Definition:** An element \\r \in C\\ is an **extreme ray** if:

- It cannot be expressed as a positive linear combination of other rays
- Exactly \\n-1\\ constraints are binding (active) at \\r\\:
  \\A\_{i_j}^T r = 0\\ for \\j \in \\1, \ldots, n-1\\\\
- The binding constraints are linearly independent

Extreme rays form the minimal generating set of the cone.

  

### Complete digraph construction

**Graph formulation:** Construct a complete directed graph \\G = (V,
E)\\ where:

- **Vertices** \\V = \\1, 2, \ldots, m\\\\ represent the \\m\\
  constraints
- **Edges** encode pairwise constraint relationships

**Associated matrix:** \\M = \[A_1, A_2, \ldots, A_m\]\\ where each
\\A_i \in \mathbb{R}^n\\ is a constraint normal vector.

  

### Spanning trees and extreme rays

**Key theorem (implicit):** There exists a correspondence between
spanning trees of the complete graph \\K_m\\ and candidate extreme rays.

**Rationale:** An extreme ray has exactly \\n-1\\ binding constraints. A
spanning tree has exactly \\n-1\\ edges (for an \\n\\-vertex tree). By
associating tree edges with binding constraints, we can systematically
enumerate extreme ray candidates.

**Complexity:**

- Number of spanning trees: \\m^{m-2}\\ (Cayley’s formula)
- Per tree: \\O(n^3)\\ for linear system solution + \\O(mn)\\ for
  feasibility check
- Total: \\O(m^{m-2} \cdot n^3)\\ — exponential in \\m\\

  

### Connection to solid angles

Once extreme rays \\\\r_1, \ldots, r_k\\\\ are obtained, the process
proceeds through simplicial decomposition (decompose the cone into
simplicial sub-cones using the extreme rays), Gram matrix computation
(for each simplicial cone, compute the Gram matrix of its normalized
extreme rays), hypergeometric series application (apply Ribando-Aomoto
formula to each simplicial cone), and inclusion-exclusion (combine solid
angles with appropriate signs).

The spanning tree method provides an **alternative to geometric
decomposition** (Brion-Vergne), grounded in combinatorial graph theory
rather than geometric algorithms.

  

## Diagnostic tools and method validation

The SolidAngleR package provides diagnostic tools to analyze cone
properties and recommend appropriate computational methods.

``` r

#### Diagnose various cones                                              ####
test_cones <- list(
  "4D orthant" = diag(4),
  "Tridiagonal" = V_trid,
  "Random 3D" = V_random
)

for (name in names(test_cones)) {
  message("--- ", name, " ---")
  diag_result <- diagnose_cone(test_cones[[name]])
  print(diag_result)
}
#> Cone diagnostic report
#> ======================
#> 
#> Dimension:                4 
#> Linearly independent:     ✓ Yes 
#> Associated matrix PD:     ✓ Yes 
#> Tridiagonal structure:    ✓ Yes 
#> 
#> Eigenvalues:              1, 1, 1, 1 
#> Minimum eigenvalue:       1 
#> 
#> Suggested method:         tridiagonal series (most efficient) ✓
#> Cone diagnostic report
#> ======================
#> 
#> Dimension:                4 
#> Linearly independent:     ✓ Yes 
#> Associated matrix PD:     ✓ Yes 
#> Tridiagonal structure:    ✓ Yes 
#> 
#> Eigenvalues:              1.114455, 1.043718, 0.956282, 0.885545 
#> Minimum eigenvalue:       0.885545 
#> 
#> Suggested method:         tridiagonal series (most efficient) ✓
#> Cone diagnostic report
#> ======================
#> 
#> Dimension:                3 
#> Linearly independent:     ✓ Yes 
#> Associated matrix PD:     ✗ No 
#> Tridiagonal structure:    ○ No 
#> 
#> Eigenvalues:              1.764742, 1.336479, -0.101221 
#> Minimum eigenvalue:       -0.101221 
#> 
#> Suggested method:         formula (closed form available) ✓
```

  

### Method comparison and validation

Different methods should yield consistent results when applicable to the
same cone.

``` r

#### 3D orthant: all methods should give 0.125                           ####
V_test <- diag(3)
expected <- 0.125

methods_3d <- c("formula", "series", "decomposition")
results_3d <- sapply(methods_3d, function(m) {
  compute_solid_angle(V_test, method = m, max_terms = 200)
})

#### Add the multivariate-normal branch of the dispatcher               ####
omega_result <- compute_solid_angle(V_test, method = "mvn")
results_3d <- c(results_3d, mvn = omega_result)

comparison_table <- data.frame(
  method = names(results_3d),
  omega_normalized = as.numeric(results_3d),
  abs_error = abs(as.numeric(results_3d) - expected)
)
knitr::kable(
  comparison_table,
  digits = c(NA, 8, 2),
  col.names = c("method", "omega (normalized)", "abs. error")
)
```

| method        | omega (normalized) | abs. error |
|:--------------|-------------------:|-----------:|
| formula       |              0.125 |          0 |
| series        |              0.125 |          0 |
| decomposition |              0.125 |          0 |
| mvn           |              0.125 |          0 |

  

## Comparing the four approaches

  

### Summary matrix

| **Method** | **Form** | **Dimensions** | **Complexity** | **Convergence** | **Exactness** |
|----|----|----|----|----|----|
| **Hypergeometric (Ribando-Aomoto)** | Infinite series | Any \\n\\ | \\O(N^{\binom{n}{2}})\\ | Positive-definite criterion | Arbitrary precision |
| **Decomposition (Fitisone-Zhou)** | Series + decomposition | Any \\n\\ | Varies; \\O(N \cdot n)\\ for tridiagonal | Guaranteed after decomposition | Arbitrary precision |
| **Geometric (Mazonka)** | Closed-form integrals | Primarily 3D | \\O(n)\\ | Not applicable (closed-form) | Exact |
| **Spanning trees (Malekmohammadi-Mostafaee)** | Graph enumeration | Any \\n\\ | \\O(m^{m-2} \cdot n^3)\\ | N/A (finds rays, not angles) | Exact rays |

  

### When to use each method

**Hypergeometric series (direct application)** should be employed when
the cone satisfies the positive-definite criterion, high-precision
requirements are present, or analytical expressions are needed; however,
this method diverges when the cone violates the criterion and becomes
impractical for very high dimensions (\\n \> 10\\).

**Decomposition methods** provide universal applicability for arbitrary
cones and exhibit major speedup when tridiagonal structure is available,
working well in moderate dimensions (\\n \leq 15\\); they incur
decomposition overhead for simple cones and face exponential blowup in
very high dimensions.

**Geometric closed-form** methods excel in three dimensions for
polyhedral cones with moderate face count, enable real-time applications
through single arctangent evaluation, and uniquely handle intersecting
cones; they do not extend to \\n \geq 4\\ and are limited to specific
geometric configurations.

**Spanning tree enumeration** converts intersection form to sum form,
systematically finds extreme rays, and supports theoretical analysis of
cone structure for small numbers of constraints (\\m \leq 10\\); it
becomes impractical for large \\m\\ due to exponential complexity and
does not compute solid angles directly.

  

### Theoretical connections

**Positive-definite criterion as unifying concept:** The criterion
determines convergence in Ribando-Aomoto theory; decomposition ensures
the criterion is met in Fitisone-Zhou methods; it is implicitly
satisfied for well-posed geometric problems in Mazonka’s approach; and
the associated matrix from found rays determines applicability of
subsequent series methods in Malekmohammadi-Mostafaee enumeration.

**Gram matrix centrality:** All methods ultimately rely on inner
products between spanning vectors; Gram matrix eigenvalues control
hypergeometric convergence rate, numerical stability, and geometric
degeneracy detection.

**Decomposition as bridge:** Decomposition connects
non-positive-definite cones to positive-definite sub-cones, enables
universal application of series methods, and parallels graph-theoretic
decomposition into spanning trees.

  

## Applications

  

### Data envelopment analysis (DEA)

**Context:** DEA measures relative efficiency of decision-making units
(DMUs) using linear programming. Weight restrictions form cones.

**Solid angle applications:** Solid angle measures the “size” of the
feasible weight space in assurance region quantification (indicating
flexibility); changes in solid angle reveal the impact of constraint
modifications in sensitivity analysis; larger solid angles for weight
restrictions imply more robust efficiency scores in efficiency dominance
assessments.

  

### Integer programming and lattice point counting

**Barvinok’s algorithm (1994)** [\[7\]](#ref7): Polynomial-time
algorithm for counting lattice points in polytopes when dimension is
fixed.

**Solid angle role:** The algorithm decomposes the polytope into tangent
cones at vertices, computes generating functions for each cone (with
solid angles appearing in coefficients), and sums contributions using
inclusion-exclusion.

**Key requirement:** Efficient solid angle computation in arbitrary
dimensions → hypergeometric methods essential.

``` r

#### Solid angle subtended by a tetrahedron face                         ####
vertices <- cbind(
  c(1, 0, 0),
  c(0, 1, 0),
  c(0, 0, 1)
)

omega_face <- solid_angle_polyhedral(vertices)
tetra_table <- data.frame(
  quantity = c("Van Oosterom formula",
               "As fraction of sphere (%)"),
  value = c(sprintf("%.6f", omega_face),
            sprintf("%.2f", omega_face * 100))
)
knitr::kable(tetra_table, caption = "Tetrahedron face solid angle.")
```

| quantity                  | value    |
|:--------------------------|:---------|
| Van Oosterom formula      | 0.125000 |
| As fraction of sphere (%) | 12.50    |

Tetrahedron face solid angle. {.table}

  

### Multivariate probability

**Orthant probabilities:** For multivariate Gaussian \\X \sim N(\mu,
\Sigma)\\:

\\P(X_1 \> 0, \ldots, X_n \> 0) =
\frac{\Omega}{(2\pi)^{n/2}\|\Sigma\|^{1/2}} + \text{correction terms}\\

Solid angle of the positive orthant under covariance-induced metric.

  

### Physics and engineering

**Radiation transport:** Solid angles determine view factors in
radiative heat transfer.

**Optics:** Light gathering power of optical systems scales with solid
angle.

**Antenna theory:** Antenna gain inversely proportional to beam solid
angle.

**Computer graphics:** Solid angles for area light sources in global
illumination.

  

### Ecological feasibility domains

A particularly important application is the quantification of
feasibility domains in multi-species ecological communities, following
Song et al. (2018) [\[8\]](#ref8).

``` r

#### 4-species community                                                 ####
S <- 4
set.seed(123)
alpha <- matrix(runif(S * S, -1, 0), S, S)
diag(alpha) <- -1

result_eco <- compute_solid_angle(alpha, method = "mvn", normalize = FALSE)

eco_table <- data.frame(
  quantity = c("Species",
               "Feasibility volume (Omega)",
               "Per-species measure",
               "Parameter-space fraction supporting coexistence (%)"),
  value = c(sprintf("%d", S),
            sprintf("%.6f", result_eco),
            sprintf("%.6f", result_eco^(1 / S)),
            sprintf("%.2f", result_eco * 100))
)
knitr::kable(eco_table, caption = "Ecological feasibility domain diagnostics (S = 4).")
```

| quantity                                            | value    |
|:----------------------------------------------------|:---------|
| Species                                             | 4        |
| Feasibility volume (Omega)                          | 0.001288 |
| Per-species measure                                 | 0.189439 |
| Parameter-space fraction supporting coexistence (%) | 0.13     |

Ecological feasibility domain diagnostics (S = 4). {.table}

  

## Conclusion

The comprehensive analysis of solid angle methods for polyhedral cones
reveals a rich mathematical landscape spanning geometry, analysis,
combinatorics, and computation. **From Girard’s 17th-century theorem**
relating angles to area, through **Euler and Lagrange’s variational
methods**, to **Van Oosterom’s efficient 3D formulation**, and finally
to **Aomoto and Ribando’s hypergeometric series** for arbitrary
dimensions, the progression demonstrates both conceptual unification and
increasing sophistication.

The four methodological approaches examined: **hypergeometric series**
[\[1\]](#ref1)[\[2\]](#ref2), **geometric closed-forms** [\[3\]](#ref3),
**decomposition algorithms** [\[2\]](#ref2), and **graph-theoretic
enumeration** [\[4\]](#ref4) provide complementary tools. Hypergeometric
methods offer analytical precision and arbitrary-dimensional
applicability but require positive-definite criteria. Fitisone and
Zhou’s decomposition breakthrough ensures universal applicability
through systematic cone subdivision, with tridiagonal optimization
providing exponential speedup. Mazonka’s geometric approach delivers
exact closed-forms and efficient complex arithmetic for 3D problems,
including the previously unsolved intersecting cones case. The spanning
tree method bridges representations, converting intersection-form cones
to sum-form via combinatorial enumeration.

**The associated matrix emerges as the unifying concept**: its
positive-definiteness determines hypergeometric convergence, its
structure guides decomposition strategies, and its relationship to graph
Laplacians connects combinatorial and analytical perspectives. The Gram
matrix of spanning vectors encodes all geometric information, with
eigenvalue spectra controlling convergence rates and numerical
stability.

Applications span data envelopment analysis (quantifying efficiency cone
flexibility), lattice point enumeration (Barvinok’s polynomial-time
algorithm), computational geometry (volume calculations, mesh quality),
multivariate statistics (orthant probabilities), and physics (radiation
view factors, antenna patterns). Modern machine learning
contexts—federated learning consensus, ensemble robustness
quantification—demonstrate continued relevance.

**Computational complexity ranges from linear** (Mazonka’s O(n)
polyhedral formula) **to exponential** (spanning trees’ \\O(m^{m-2})\\,
general hypergeometric’s \\O(N^{\binom{n}{2}})\\), with tridiagonal
structure providing the crucial middle ground enabling practical
high-dimensional computation.

This synthesis demonstrates that **solid angle theory exemplifies
mathematical maturity**: classical results guide intuition, modern
analysis ensures rigor, computational methods enable application, and
ongoing research pushes boundaries. The interplay between geometry and
analysis, continuous and discrete mathematics, theory and computation,
makes solid angle methods an enduring subject of mathematical
investigation with practical impact across scientific disciplines.

  

## References

**\[1\]** Ribando, J. M. (2006). Measuring solid angles beyond dimension
three. *Discrete & Computational Geometry*, 36(3), 479-487. DOI:
[10.1007/s00454-006-1253-4](https://doi.org/10.1007/s00454-006-1253-4)

**\[2\]** Fitisone, A., & Zhou, Y. (2023). Solid angle measure of
polyhedral cones. arXiv:2304.11102 \[math.CO\]. URL:
<https://arxiv.org/abs/2304.11102>

**\[3\]** Mazonka, O. (2012). Solid angle of conical surfaces,
polyhedral cones, and intersecting spherical caps. arXiv:1205.1396v2
\[math.MG\]. URL: <https://arxiv.org/abs/1205.1396>

**\[4\]** Malekmohammadi, N., & Mostafaee, A. (2017). Obtaining all
extreme rays of a special cone using spanning trees in a complete
digraph: application in DEA. *Journal of the Operational Research
Society*, 69(3), 465-472. DOI:
[10.1057/s41274-017-0265-9](https://doi.org/10.1057/s41274-017-0265-9)

**\[5\]** Van Oosterom, A., & Strackee, J. (1983). The solid angle of a
plane triangle. *IEEE Transactions on Biomedical Engineering*,
BME-30(2), 125-126. DOI:
[10.1109/TBME.1983.325207](https://doi.org/10.1109/TBME.1983.325207)

**\[6\]** Aomoto, K. (1977). Analytic structure of Schläfli function.
*Nagoya Mathematical Journal*, 68, 1-16. URL:
<https://projecteuclid.org/euclid.nmj/1118786429>

**\[7\]** Barvinok, A. I. (1994). A polynomial time algorithm for
counting integral points in polyhedra when the dimension is fixed.
*Mathematics of Operations Research*, 19(4), 769-779. DOI:
[10.1287/moor.19.4.769](https://doi.org/10.1287/moor.19.4.769)

**\[8\]** Song, C., Rohr, R. P., & Saavedra, S. (2018). A guideline to
study the feasibility domain of multi-trophic and changing ecological
communities. *Journal of Theoretical Biology*, 450, 30-36. DOI:
[10.1016/j.jtbi.2018.04.030](https://doi.org/10.1016/j.jtbi.2018.04.030)

**\[9\]** Girard, A. (1626). *Tables des sinus, tangentes, & secantes.*
Leiden.

**\[10\]** Euler, L. (1755). Principes de la trigonométrie sphérique
tirés de la méthode des plus grands et des plus petits. *Mémoires de
l’Académie des Sciences de Berlin*, 9, 223-257.

**\[11\]** L’Huilier, S. A. J. (1789). *Exposition élémentaire des
principes des calculs supérieurs.* Berlin.

**\[12\]** Schläfli, L. (1858-1860). On the multiple integral…
*Quarterly Journal of Pure and Applied Mathematics*, 3, 54-68, 97-108.
