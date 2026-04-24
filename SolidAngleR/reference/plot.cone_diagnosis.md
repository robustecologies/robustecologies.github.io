# Plot the eigenvalue spectrum of a diagnosed cone

Visualises the eigenvalue spectrum of the associated matrix \\M_n(C)\\
as a signed bar chart. Positive eigenvalues are coloured blue; negative
eigenvalues, which signal that the dispatcher must fall back to the
decomposition method, are coloured red. A dashed reference line is drawn
at zero.

## Usage

``` r
# S3 method for class 'cone_diagnosis'
plot(x, type = "eigenvalues", ...)
```

## Arguments

- x:

  An object of class `cone_diagnosis` returned by
  [`diagnose_cone`](https://robustecologies.github.io/SolidAngleR/reference/diagnose_cone.md).

- type:

  Character. Currently only `"eigenvalues"` is supported; future plot
  types will be added without breaking this argument.

- ...:

  Further arguments passed to
  [`barplot`](https://rdrr.io/r/graphics/barplot.html).

## Value

`NULL`, invisibly. Called for the side effect of producing the plot.

## Details

The eigenvalues of \\M_n(C)\\ are the discriminant of method selection
(Ribando 2006, Theorem 1.5; Fitisone & Zhou 2023, Theorem 3.3). When all
eigenvalues are positive, the hypergeometric series converges and is the
preferred backend; the appearance of any negative eigenvalue forces the
decomposition route. The plot makes that diagnostic visible at a glance.

## References

Ribando, J. M. (2006). Measuring solid angles beyond dimension three.
*Discrete & Computational Geometry*, 36(3), 479-487.
[doi:10.1007/s00454-006-1253-4](https://doi.org/10.1007/s00454-006-1253-4)

Fitisone, A., & Zhou, Y. (2023). Solid angle measure of polyhedral
cones. arXiv:2304.11102 (math.CO).

## See also

[`diagnose_cone`](https://robustecologies.github.io/SolidAngleR/reference/diagnose_cone.md)
for the constructor;
[`print.cone_diagnosis`](https://robustecologies.github.io/SolidAngleR/reference/print.cone_diagnosis.md),
[`summary.cone_diagnosis`](https://robustecologies.github.io/SolidAngleR/reference/summary.cone_diagnosis.md)
for the companion S3 methods;
[`compute_solid_angle`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angle.md)
for the dispatcher whose method choice this plot informs.

## Examples

``` r
if (FALSE) { # \dontrun{
# Positive-definite cone (orthant)
V <- diag(4)
plot(diagnose_cone(V))

# Non-positive-definite cone
V <- matrix(c(1, 0.9, 0.1,
              0, 0.5, 0.8,
              0, 0, 0.3), nrow = 3)
plot(diagnose_cone(V))
} # }
```
