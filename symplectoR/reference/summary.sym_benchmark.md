# Summarize a benchmark problem

The print view plus generator metadata (condition number, damping
constants, provenance).

## Usage

``` r
# S3 method for class 'sym_benchmark'
summary(object, ...)

# S3 method for class 'summary.sym_benchmark'
print(x, ...)
```

## Arguments

- object:

  A `sym_benchmark` object.

- ...:

  Ignored.

- x:

  A `summary.sym_benchmark` object.

## Value

An object of class `summary.sym_benchmark`.

## See also

[`sym_benchmark()`](https://robustecologies.github.io/symplectoR/reference/sym_benchmark.md),
[`print.sym_benchmark()`](https://robustecologies.github.io/symplectoR/reference/print.sym_benchmark.md),
[`plot.sym_benchmark()`](https://robustecologies.github.io/symplectoR/reference/plot.sym_benchmark.md)

## Examples

``` r
if (FALSE) { # \dontrun{
summary(sym_benchmark("damped_oscillator"))
} # }
```
