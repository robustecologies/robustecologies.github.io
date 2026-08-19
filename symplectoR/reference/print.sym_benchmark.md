# Print a benchmark problem

One-screen report of a
[`sym_benchmark()`](https://robustecologies.github.io/symplectoR/reference/sym_benchmark.md):
dimension, box, ground-truth minimum and analytic assets.

## Usage

``` r
# S3 method for class 'sym_benchmark'
print(x, ...)
```

## Arguments

- x:

  A `sym_benchmark` object.

- ...:

  Ignored.

## Value

`x`, invisibly.

## See also

[`sym_benchmark()`](https://robustecologies.github.io/symplectoR/reference/sym_benchmark.md),
[`summary.sym_benchmark()`](https://robustecologies.github.io/symplectoR/reference/summary.sym_benchmark.md),
[`plot.sym_benchmark()`](https://robustecologies.github.io/symplectoR/reference/plot.sym_benchmark.md)

## Examples

``` r
if (FALSE) { # \dontrun{
print(sym_benchmark("quadratic", d = 10, kappa = 100))
} # }
```
