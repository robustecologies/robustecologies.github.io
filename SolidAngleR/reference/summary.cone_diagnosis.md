# Extended diagnostic summary for a polyhedral cone

Computes spectral diagnostics beyond the basic `print`: condition number
of the associated matrix, eigenvalue spread, a five-level stability
classification, and a method-selection rationale that records the
discriminating quantity (typically \\\kappa(M_n(C))\\) used by the
dispatcher.

## Usage

``` r
# S3 method for class 'cone_diagnosis'
summary(object, ...)
```

## Arguments

- object:

  An object of class `cone_diagnosis` returned by
  [`diagnose_cone`](https://robustecologies.github.io/SolidAngleR/reference/diagnose_cone.md).

- ...:

  Currently unused; present for S3 generic compatibility.

## Value

An object of class `summary.cone_diagnosis`; a list with components

- dimension:

  Ambient dimension of the cone.

- eigenvalues:

  Eigenvalues of the associated matrix.

- condition_number:

  Ratio of the largest to the smallest absolute eigenvalue, or `Inf`
  when the smallest is below `.Machine$double.eps`.

- eigenvalue_spread:

  Difference between the largest and smallest eigenvalue.

- stability:

  One of `"excellent"`, `"good"`, `"moderate"`, `"poor"`, `"critical"`,
  derived from the condition number with thresholds at \\10\\, \\10^2\\,
  \\10^6\\, \\10^{10}\\.

- method_rationale:

  Character string explaining the recommended backend in terms of the
  spectrum.

- original:

  The original `cone_diagnosis` object.

## Details

The condition number is computed on the absolute spectrum to remain
meaningful when \\M_n(C)\\ is not positive semidefinite (in which case
the dispatcher uses the decomposition route). The five stability bands
are intentionally coarse and serve only to flag situations where the
user should expect numerical sensitivity in the chosen backend.

## References

Ribando, J. M. (2006). Measuring solid angles beyond dimension three.
*Discrete & Computational Geometry*, 36(3), 479-487.
[doi:10.1007/s00454-006-1253-4](https://doi.org/10.1007/s00454-006-1253-4)

## See also

[`diagnose_cone`](https://robustecologies.github.io/SolidAngleR/reference/diagnose_cone.md)
for the constructor;
[`print.cone_diagnosis`](https://robustecologies.github.io/SolidAngleR/reference/print.cone_diagnosis.md),
[`plot.cone_diagnosis`](https://robustecologies.github.io/SolidAngleR/reference/plot.cone_diagnosis.md)
for the companion S3 methods;
[`print.summary.cone_diagnosis`](https://robustecologies.github.io/SolidAngleR/reference/print.summary.cone_diagnosis.md)
for the print method on the summary object.

## Examples

``` r
if (FALSE) { # \dontrun{
V <- diag(4)
d <- diagnose_cone(V)
s <- summary(d)
print(s)
} # }
```
