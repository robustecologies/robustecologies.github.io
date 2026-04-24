# Print a cone diagnosis report

Prints a one-screen synopsis of a `cone_diagnosis` object: ambient
dimension, linear independence flag, positive definiteness of the
associated matrix, tridiagonal structure, eigenvalue spectrum, and the
recommended computation method.

## Usage

``` r
# S3 method for class 'cone_diagnosis'
print(x, ...)
```

## Arguments

- x:

  An object of class `cone_diagnosis` returned by
  [`diagnose_cone`](https://robustecologies.github.io/SolidAngleR/reference/diagnose_cone.md).

- ...:

  Currently unused; present for S3 generic compatibility.

## Value

The input object `x` returned invisibly.

## Details

Side effect only: writes a formatted report to the console. Use
[`summary.cone_diagnosis`](https://robustecologies.github.io/SolidAngleR/reference/summary.cone_diagnosis.md)
for the extended report including condition number, eigenvalue spread,
stability classification, and the rationale behind the method
recommendation.

## See also

[`diagnose_cone`](https://robustecologies.github.io/SolidAngleR/reference/diagnose_cone.md)
for the constructor;
[`summary.cone_diagnosis`](https://robustecologies.github.io/SolidAngleR/reference/summary.cone_diagnosis.md),
[`plot.cone_diagnosis`](https://robustecologies.github.io/SolidAngleR/reference/plot.cone_diagnosis.md)
for the companion S3 methods.

## Examples

``` r
if (FALSE) { # \dontrun{
V <- diag(4)
d <- diagnose_cone(V)
print(d)
} # }
```
