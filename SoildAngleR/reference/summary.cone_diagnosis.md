# Summary method for cone diagnosis

Provides extended diagnostic information for a polyhedral cone,
including spectral analysis of the associated matrix, condition number
assessment, and detailed method selection rationale.

## Usage

``` r
# S3 method for class 'cone_diagnosis'
summary(object, ...)
```

## Arguments

- object:

  A cone_diagnosis object from
  [`diagnose_cone`](https://robustecologies.github.io/SolidAngleR/reference/diagnose_cone.md)

- ...:

  Additional arguments (unused)

## Value

An object of class `summary.cone_diagnosis` containing:

- dimension:

  Ambient dimension of the cone

- eigenvalues:

  Eigenvalues of the associated matrix

- condition_number:

  Condition number (ratio of max to min eigenvalue)

- eigenvalue_spread:

  Difference between maximum and minimum eigenvalues

- stability:

  Stability assessment based on condition number

- method_rationale:

  Detailed explanation for method recommendation

- original:

  The original cone_diagnosis object

## Examples

``` r
if (FALSE) { # \dontrun{
# Create a test cone
V <- diag(4)

# Diagnose and summarize
d <- diagnose_cone(V)
print(d)
s <- summary(d)
print(s)

# Visualize eigenvalue structure
plot(d)
} # }
```
