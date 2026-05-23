# Print method for a summary.solver_spec object

Formats the solver-taxonomy summary as a compact report: method family,
consistency order, stability characterisation and the principal control
parameters.

## Usage

``` r
# S3 method for class 'summary.solver_spec'
print(x, ...)
```

## Arguments

- x:

  A `summary.solver_spec` object.

- ...:

  Unused.

## Value

The input `x`, invisibly.

## See also

[`summary.solver_spec()`](https://robustecologies.github.io/janos/reference/summary.solver_spec.md)
. summary constructor.

## Examples

``` r
if (FALSE) { # \dontrun{
print(summary(solver_rk45()))
} # }
```
