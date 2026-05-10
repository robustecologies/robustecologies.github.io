# Logical check of data

This function provides a logical test of whether or not a given list of
data meets minimum criteria to be considered data for
[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md), or
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Usage

``` r
is.data(Data)
```

## Arguments

- Data:

  A list of data. For more information, see the
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
  function.

## Value

A logical value. Returns `TRUE` if `Data` meets minimum requirements to
be considered data in this package, and `FALSE` otherwise.

## Details

This function is useful for testing whether or not a list of data meets
minimum criteria to be considered data in this package. The minimum
requirements are that `Data` is a list, and it contains `mon.names` and
`parm.names`.

This function is not extensive. For example, it does not match the
length of `parm.names` with the length of `Initial.Values`, or compare
the length of `mon.names` to the number of monitored variables output
from the `Model` specification function. Additional checks are conducted
in
[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md), and
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## See also

[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Valid data list
Data <- list(N = 25, y = rnorm(25), mon.names = "LP",
  parm.names = c("mu", "log.sigma"))
is.data(Data) # TRUE

# Missing required components
Bad <- list(y = rnorm(25))
is.data(Bad) # FALSE
} # }
```
