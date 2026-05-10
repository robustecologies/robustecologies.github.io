# Prior predictive check

Generates prior predictive data by sampling from the prior distribution
and evaluating the model, then compares the predictive distribution to
the observed data.

## Usage

``` r
prior_predictive_check(
  Model,
  Data,
  rprior,
  n = 500,
  type = "density",
  stat_fun = mean,
  col = NULL,
  ...
)
```

## Arguments

- Model:

  The model specification function, as supplied to
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

- Data:

  The data list, as supplied to
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

- rprior:

  A function that takes a single integer `n` and returns an \\n \times
  K\\ matrix of prior draws (one row per sample, one column per
  parameter).

- n:

  Integer; number of prior samples to draw (default 500).

- type:

  Character string selecting the plot type. One of `"density"`
  (default), `"intervals"`, `"ribbon"`, or `"stat"`.

- stat_fun:

  A function computing a scalar test statistic from a numeric vector.
  Used when `type = "stat"` (default: `mean`).

- col:

  Optional character vector of colors. When `NULL` (default), the RElab
  palette is used.

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns the ggplot object.

## Details

The function evaluates the model at each prior draw to obtain predicted
values, then overlays the prior predictive distribution on the observed
data. This helps diagnose whether the prior is unreasonably diffuse or
informative relative to the data.

## References

Gabry, J., Simpson, D., Vehtari, A., Betancourt, M., and Gelman, A.
(2019). "Visualization in Bayesian workflow." *Journal of the Royal
Statistical Society: Series A*, 182(2), 389–402.
[doi:10.1111/rssa.12378](https://doi.org/10.1111/rssa.12378)

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
rprior <- function(n) {
    cbind(rnorm(n, 0, 10), rexp(n, 1))
}
prior_predictive_check(Model, Data, rprior, n = 200)
} # }
```
