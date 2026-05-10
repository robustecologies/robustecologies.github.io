# Plot method for SMC objects

Produces diagnostic plots for an SMC fit. The default type shows
weighted marginal posterior histograms. Alternative types include ESS
trajectory, tempering schedule, intervals, and pairs.

## Usage

``` r
# S3 method for class 'smc'
plot(x, type = "posterior", ground_truth = NULL, col = NULL, ...)
```

## Arguments

- x:

  An object of class `smc`.

- type:

  Character: `"posterior"` (default), `"ess"`, `"schedule"`,
  `"intervals"`, or `"pairs"`.

- ground_truth:

  Optional named numeric vector of true parameter values. When provided
  and `type = "posterior"`, vertical dashed lines are drawn at the true
  values.

- col:

  Optional character vector of hex color strings. When non-`NULL`,
  overrides the default RElab contrasting palette.

- ...:

  Additional arguments passed to plotting functions.

## Value

Invisibly returns the plot(s).

## Details

Produces diagnostic plots of a sequential Monte Carlo fit produced by
[`SMC`](https://robustecologies.github.io/lucifer/reference/SMC.md).
Summary of the content is given below. Default output renders a
multi-panel graphic (trace, density, and autocorrelation where
applicable). The `PDF` argument captures the graphic to a file;
otherwise the current device is used. Font and colour choices follow
[`theme_relab`](https://robustecologies.github.io/lucifer/reference/theme_relab.md).

## References

Del Moral, P., Doucet, A., & Jasra, A. (2006). Sequential Monte Carlo
samplers. *Journal of the Royal Statistical Society, Series B*, 68(3),
411-436.
[doi:10.1111/j.1467-9868.2006.00553.x](https://doi.org/10.1111/j.1467-9868.2006.00553.x)

## See also

[`print.smc`](https://robustecologies.github.io/lucifer/reference/print.smc.md),
[`summary.smc`](https://robustecologies.github.io/lucifer/reference/summary.smc.md),
[`to_mcmc_list.smc`](https://robustecologies.github.io/lucifer/reference/to_mcmc_list.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving plot.smc
} # }
```
