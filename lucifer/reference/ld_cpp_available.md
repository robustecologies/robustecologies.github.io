# Check C++ backend availability

Returns the number of OpenMP threads available to the C++ backend. A
return value of 1 indicates either single-threaded operation or that the
package was compiled without OpenMP support.

## Usage

``` r
ld_cpp_available()
```

## Value

Integer number of available OpenMP threads.

## Details

Implementation of `ld_cpp_available`. Refer to the package vignettes and
the cited references for a complete algorithmic and mathematical
description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
ld_cpp_available()
} # }
```
