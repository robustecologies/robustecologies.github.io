# Print an objective specification

One-screen report of a
[`sym_objective()`](https://robustecologies.github.io/symplectoR/reference/sym_objective.md):
dimension, evaluation route, gradient availability and constraint
handling.

## Usage

``` r
# S3 method for class 'sym_objective'
print(x, ...)
```

## Arguments

- x:

  A `sym_objective` object.

- ...:

  Ignored.

## Value

`x`, invisibly.

## See also

[`sym_objective()`](https://robustecologies.github.io/symplectoR/reference/sym_objective.md),
[`summary.sym_objective()`](https://robustecologies.github.io/symplectoR/reference/summary.sym_objective.md),
[`plot.sym_objective()`](https://robustecologies.github.io/symplectoR/reference/plot.sym_objective.md)

## Examples

``` r
if (FALSE) { # \dontrun{
print(sym_objective(function(x) sum(x^2), dim = 3))
} # }
```
