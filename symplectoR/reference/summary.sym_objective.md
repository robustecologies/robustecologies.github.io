# Summarize an objective specification

Reports the print view plus a probe evaluation at the box midpoint (or
the origin), confirming the objective is callable.

## Usage

``` r
# S3 method for class 'sym_objective'
summary(object, ...)

# S3 method for class 'summary.sym_objective'
print(x, ...)
```

## Arguments

- object:

  A `sym_objective` object.

- ...:

  Ignored.

- x:

  A `summary.sym_objective` object.

## Value

An object of class `summary.sym_objective`.

## See also

[`sym_objective()`](https://robustecologies.github.io/symplectoR/reference/sym_objective.md),
[`print.sym_objective()`](https://robustecologies.github.io/symplectoR/reference/print.sym_objective.md),
[`plot.sym_objective()`](https://robustecologies.github.io/symplectoR/reference/plot.sym_objective.md)

## Examples

``` r
if (FALSE) { # \dontrun{
summary(sym_objective(function(x) sum(x^2), dim = 3, lower = -1, upper = 1))
} # }
```
