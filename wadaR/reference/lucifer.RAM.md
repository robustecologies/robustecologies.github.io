# lucifer RAM estimate

Estimates the random-access memory (RAM) required to update a given
model and data with the
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
function.

## Usage

``` r
lucifer.RAM(Model, Data, Iterations, Thinning, Algorithm = "RWM")
```

## Arguments

- Model:

  A model specification function. See
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

- Data:

  A list of data. See
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

- Iterations:

  The number of iterations for which
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
  would update.

- Thinning:

  The amount of thinning applied to the chains.

- Algorithm:

  The name of the algorithm as a string, as entered in
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).
  Defaults to `"RWM"`.

## Value

A list with estimated RAM (in MB) for each component: Covar, Data,
Deviance, Initial.Values, Model, Monitor, post, Posterior1, Posterior2,
Summary1, Summary2, and Total.

## Details

The `lucifer.RAM` function uses
[`object.size`](https://rdrr.io/r/utils/object.size.html) to estimate
the size in MB of RAM required to update one chain in
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
for a given model and data. When RAM is exceeded, the computer will
crash. This function is useful when estimating how many iterations to
update without crashing.

## See also

[`BigData`](https://robustecologies.github.io/lucifer/reference/BigData.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LML`](https://robustecologies.github.io/lucifer/reference/LML.md),
[`object.size`](https://rdrr.io/r/utils/object.size.html)

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving lucifer.RAM
} # }
```
