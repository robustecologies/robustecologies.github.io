# Combine demonoid objects

This function combines objects of class `demonoid`.

## Usage

``` r
Combine(x, Data, Thinning = 1, BurnIn = NULL)
```

## Arguments

- x:

  A list of objects of class `demonoid`.

- Data:

  The data list, which must be identical to the data used to create the
  `demonoid` objects with
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

- Thinning:

  The amount of thinning to apply to the posterior samples after
  appending them together. Defaults to 1.

- BurnIn:

  Optional integer specifying the number of warm-up iterations to
  discard from each chain before combining. If `NULL` (the default), no
  extra burn-in is applied beyond what each chain already carries.

## Value

An object of class `demonoid`. For more information, see
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

## Details

The `Combine` function enables a user to combine objects of class
`demonoid` for one of three reasons. Parallel chains from
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
(with `Chains > 1`) may be combined after convergence is assessed with
[`Gelman.Diagnostic`](https://robustecologies.github.io/lucifer/reference/Gelman.Diagnostic.md).
Consecutive updates of single chains may be combined when the computer
has insufficient RAM for the user to update once with enough iterations.
Multiple updates may be combined when LP oscillates, suggesting that
thinning is insufficient.

## See also

[`caterpillar.plot`](https://robustecologies.github.io/lucifer/reference/caterpillar.plot.md),
[`Gelman.Diagnostic`](https://robustecologies.github.io/lucifer/reference/Gelman.Diagnostic.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`Thin`](https://robustecologies.github.io/lucifer/reference/Thin.md)

## Examples

``` r
if (FALSE) { # \dontrun{
Fit3 <- Combine(list(Fit1, Fit2), Data=MyData, Thinning=2)
} # }
```
