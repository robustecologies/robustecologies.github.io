# Display parallelization diagnostic information

Prints a diagnostic report of the current parallelization configuration
including detected hardware, OpenMP availability, current thread budget,
and platform-specific guidance.

## Usage

``` r
lucifer_parallel_info()
```

## Value

Invisibly returns the diagnostic information as a list.

## Details

Implementation of `lucifer_parallel_info`. Refer to the package
vignettes and the cited references for a complete algorithmic and
mathematical description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
lucifer_parallel_info()
} # }
```
