# Print a parameter sweep

One-screen report of a
[`sym_sweep()`](https://robustecologies.github.io/symplectoR/reference/sym_sweep.md)
run: grid size, success fraction and the incumbent cell.

## Usage

``` r
# S3 method for class 'sym_sweep'
print(x, ...)
```

## Arguments

- x:

  A `sym_sweep` object.

- ...:

  Ignored.

## Value

`x`, invisibly.

## See also

[`sym_sweep()`](https://robustecologies.github.io/symplectoR/reference/sym_sweep.md),
[`summary.sym_sweep()`](https://robustecologies.github.io/symplectoR/reference/summary.sym_sweep.md),
[`plot.sym_sweep()`](https://robustecologies.github.io/symplectoR/reference/plot.sym_sweep.md)

## Examples

``` r
if (FALSE) { # \dontrun{
sw <- sym_sweep(sym_objective(sym_benchmark("quadratic", d = 4)), "slc_expo",
                grid = list(C = 10^(-2:2), h = 10^(-2:1)), x0 = rep(2, 4))
print(sw)
} # }
```
