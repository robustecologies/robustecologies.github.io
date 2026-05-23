# As posterior predictive check

This function converts an object of class `demonoid.val` to an object of
class `demonoid.ppc`.

## Usage

``` r
as.ppc(x, set = 3)
```

## Arguments

- x:

  This is an object of class `demonoid.val`.

- set:

  This is an integer that indicates which list component is to be used.
  When `set=1`, the modeled data set is used. When `set=2`, the
  validation data set is used. When `set=3`, both data sets are used.

## Value

The returned object is an object of class `demonoid.ppc`.

## Details

After using the
[`Validate`](https://robustecologies.github.io/lucifer/reference/Validate.md)
function for holdout validation, it is often suggested to perform
posterior predictive checks. The `as.ppc` function converts the output
object of
[`Validate`](https://robustecologies.github.io/lucifer/reference/Validate.md),
which is an object of class `demonoid.val`, to an object of class
`demonoid.ppc`. The returned object is the same as if it were created
with the
[`predict.demonoid`](https://robustecologies.github.io/lucifer/reference/predict.demonoid.md)
function, rather than the
[`Validate`](https://robustecologies.github.io/lucifer/reference/Validate.md)
function.

After this conversion, the user may use posterior predictive checks, as
usual, with the
[`summary.demonoid.ppc`](https://robustecologies.github.io/lucifer/reference/summary.demonoid.ppc.md)
function.

## See also

[`predict.demonoid`](https://robustecologies.github.io/lucifer/reference/predict.demonoid.md),
[`summary.demonoid.ppc`](https://robustecologies.github.io/lucifer/reference/summary.demonoid.ppc.md),
and
[`Validate`](https://robustecologies.github.io/lucifer/reference/Validate.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving as.ppc
} # }
```
