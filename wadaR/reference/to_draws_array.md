# Convert lucifer fits to posterior draws formats

Converts a `demonoid` object to the `draws_array`, `draws_matrix`, or
`draws_df` format from the posterior package.

S3 method: apply `to_draws_array()` to objects of class `demonoid`.

See Details.

S3 method: apply `to_draws_matrix()` to objects of class `demonoid`.

See Details.

S3 method: apply `to_draws_df()` to objects of class `demonoid`.

## Usage

``` r
to_draws_array(x, ...)

# S3 method for class 'demonoid'
to_draws_array(x, ...)

to_draws_matrix(x, ...)

# S3 method for class 'demonoid'
to_draws_matrix(x, ...)

to_draws_df(x, ...)

# S3 method for class 'demonoid'
to_draws_df(x, ...)
```

## Arguments

- x:

  A `demonoid` object.

- ...:

  Additional arguments (currently unused).

## Value

An object of the corresponding posterior draws class.

See Details.

See Details.

See Details.

See Details.

See Details.

## Details

Implementation of `to_draws_array`. Refer to the package vignettes and
the cited references for a complete algorithmic and mathematical
description.

Implementation of `to_draws_array.demonoid`. Refer to the package
vignettes and the cited references for a complete algorithmic and
mathematical description.

Implementation of `to_draws_matrix`. Refer to the package vignettes and
the cited references for a complete algorithmic and mathematical
description.

Implementation of `to_draws_matrix.demonoid`. Refer to the package
vignettes and the cited references for a complete algorithmic and
mathematical description.

Implementation of `to_draws_df`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `to_draws_df.demonoid`. Refer to the package vignettes
and the cited references for a complete algorithmic and mathematical
description.

## See also

[`as.demonoid`](https://robustecologies.github.io/lucifer/reference/as.demonoid.md),
[`to_mcmc_list`](https://robustecologies.github.io/lucifer/reference/to_mcmc_list.md)

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LOO.demonoid`](https://robustecologies.github.io/lucifer/reference/LOO.md),
[`WAIC.demonoid`](https://robustecologies.github.io/lucifer/reference/WAIC.md),
[`log_lik.demonoid`](https://robustecologies.github.io/lucifer/reference/log_lik.md),
[`plot.demonoid`](https://robustecologies.github.io/lucifer/reference/plot.demonoid.md),
[`predict.demonoid`](https://robustecologies.github.io/lucifer/reference/predict.demonoid.md),
[`print.demonoid`](https://robustecologies.github.io/lucifer/reference/print.demonoid.md).

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- lucifer(Model, Data, IV, Algorithm = "NUTS")
draws <- to_draws_array(fit)
posterior::summarise_draws(draws)
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving to_draws_array.demonoid
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving to_draws_matrix
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving to_draws_matrix.demonoid
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving to_draws_df
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving to_draws_df.demonoid
} # }
```
