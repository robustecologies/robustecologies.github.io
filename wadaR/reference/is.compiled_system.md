# Check if an object is a compiled_system

Predicate that tests S3 inheritance against the `compiled_system` class.
Returns `TRUE` for objects produced by
[`compiled_system()`](https://robustecologies.github.io/wadaR/reference/compiled_system.md)
and `FALSE` for any other input, including `NULL`, base R objects, and
the legacy `dynamical_system` class that is no longer supported.

## Usage

``` r
is.compiled_system(x)
```

## Arguments

- x:

  Object to test.

## Value

Logical scalar. `TRUE` if `x` inherits from `compiled_system` and
therefore carries a pre-compiled C++/OpenMP basin function in
`$compiled_basin_func`; `FALSE` otherwise.

## Details

Use this predicate in user code that branches on whether a system has
been processed through the C++ compilation pipeline. Internally,
[`compute_basins()`](https://robustecologies.github.io/wadaR/reference/compute_basins.md)
dispatches to a different code path for compiled systems versus legacy
R-only system objects, and downstream functions such as
[`bifurcation_basins()`](https://robustecologies.github.io/wadaR/reference/bifurcation_basins.md)
require a compiled system.

## See also

[`compiled_system()`](https://robustecologies.github.io/wadaR/reference/compiled_system.md)
for the constructor;
[`print.compiled_system()`](https://robustecologies.github.io/wadaR/reference/print.compiled_system.md),
[`summary.compiled_system()`](https://robustecologies.github.io/wadaR/reference/summary.compiled_system.md)
and
[`plot.compiled_system()`](https://robustecologies.github.io/wadaR/reference/plot.compiled_system.md)
for the S3 method trio;
[`compile_basin_function()`](https://robustecologies.github.io/wadaR/reference/compile_basin_function.md)
for the lower-level compilation entry point.

## Examples

``` r
if (FALSE) { # \dontrun{
duffing <- compiled_system(
    cpp_dynamics = "deriv[0] = state[1]; deriv[1] = -0.3*state[1] - state[0];",
    attractors = list(
        attractor_point(c(1, 0), 0.3),
        attractor_point(c(-1, 0), 0.3),
        attractor_point(c(0, 0), 0.3)
    )
)
is.compiled_system(duffing)   # TRUE
is.compiled_system(list())    # FALSE
} # }
```
