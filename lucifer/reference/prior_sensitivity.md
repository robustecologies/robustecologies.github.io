# Prior sensitivity visualization

Overlays posterior densities obtained under different prior
specifications to visualize prior sensitivity.

## Usage

``` r
prior_sensitivity(fits, Parms = NULL, ground_truth = NULL, col = NULL, ...)
```

## Arguments

- fits:

  A named list of fitted objects (all of the same class). Names are used
  as legend labels.

- Parms:

  Character vector of parameter names to select (uses `grep` matching).
  Defaults to all parameters (up to 9).

- ground_truth:

  Optional named numeric vector of true parameter values.

- col:

  Optional character vector of colors. When `NULL` (default), the RElab
  palette is used.

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns the list of ggplot objects.

## Details

Implementation of `prior_sensitivity`. Refer to the package vignettes
and the cited references for a complete algorithmic and mathematical
description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
prior_sensitivity(
    list("N(0,1)"  = fit_narrow,
         "N(0,10)" = fit_wide,
         "N(0,100)" = fit_vague),
    ground_truth = c("beta[1]" = 2)
)
} # }
```
