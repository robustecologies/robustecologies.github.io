# Print an extended cone diagnostic summary

Prints the extended report produced by
[`summary.cone_diagnosis`](https://robustecologies.github.io/SolidAngleR/reference/summary.cone_diagnosis.md):
basic properties, spectral analysis (eigenvalues, condition number,
spread, stability band) and the method-selection rationale.

## Usage

``` r
# S3 method for class 'summary.cone_diagnosis'
print(x, ...)
```

## Arguments

- x:

  An object of class `summary.cone_diagnosis`.

- ...:

  Currently unused; present for S3 generic compatibility.

## Value

The input object `x` returned invisibly.

## Details

Side effect only: writes the formatted report to the console using the
three-section layout (basic properties, spectral analysis, method
rationale).

## See also

[`summary.cone_diagnosis`](https://robustecologies.github.io/SolidAngleR/reference/summary.cone_diagnosis.md)
for the constructor of the summary object;
[`diagnose_cone`](https://robustecologies.github.io/SolidAngleR/reference/diagnose_cone.md),
[`print.cone_diagnosis`](https://robustecologies.github.io/SolidAngleR/reference/print.cone_diagnosis.md),
[`plot.cone_diagnosis`](https://robustecologies.github.io/SolidAngleR/reference/plot.cone_diagnosis.md)
for the rest of the S3 family.

## Examples

``` r
if (FALSE) { # \dontrun{
V <- diag(4)
d <- diagnose_cone(V)
s <- summary(d)
print(s)
} # }
```
