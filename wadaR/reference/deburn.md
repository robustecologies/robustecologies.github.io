# De-burn

The `deburn` function discards or removes a user-specified number of
burn-in iterations from an object of class `demonoid`.

## Usage

``` r
deburn(x, BurnIn = 0)
```

## Arguments

- x:

  This is an object of class `demonoid`.

- BurnIn:

  This argument defaults to `BurnIn=0`, and accepts an integer that
  indicates the number of iterations to discard as burn-in.

## Value

The `deburn` function returns an object of class `demonoid`.

## Details

Documentation for the
[`burnin`](https://robustecologies.github.io/lucifer/reference/burnin.md)
function provides an introduction to the concept of burn-in as it
relates to Markov chains.

The `deburn` function discards a number of the first posterior samples,
as specified by the `BurnIn` argument. Stationarity is not checked,
because it is assumed the user has a reason for using the `deburn`
function, rather than using the results from the object of class
`demonoid`. Therefore, the posterior samples in `Posterior1` and
`Posterior2` are identical, as are `Summary1` and `Summary2`.

## See also

[`burnin`](https://robustecologies.github.io/lucifer/reference/burnin.md)
and
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

## Examples

``` r
if (FALSE) { # \dontrun{
### Assuming the user has Fit which is an object of class demonoid:
Fit2 <- deburn(Fit, BurnIn=100)
} # }
```
