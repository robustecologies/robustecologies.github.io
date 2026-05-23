# Plot method for SSM fit objects

Visualizes state-space model inference results.

## Usage

``` r
# S3 method for class 'ssm_fit'
plot(x, type = c("states", "parameters", "pairs"), state = 1L, obs = TRUE, ...)
```

## Arguments

- x:

  An object of class `ssm_fit`.

- type:

  Character string specifying the plot type: `"states"` (default, state
  trajectories with credible bands), `"parameters"` (trace plots and
  densities for static parameters), `"pairs"` (bivariate scatterplots).

- state:

  Integer index of the state variable to plot when `type = "states"`
  (default 1).

- obs:

  Logical; if `TRUE` (default), overlay observed data as points on state
  trajectory plots. Only used when `type = "states"`.

- ...:

  Additional arguments passed to plotting functions.

## Value

Invisibly returns `NULL`.

## Details

Produces diagnostic plots of a state-space model fit produced by
[`SSM`](https://robustecologies.github.io/lucifer/reference/SSM.md).
Summary of the content is given below. Default output renders a
multi-panel graphic (trace, density, and autocorrelation where
applicable). The `PDF` argument captures the graphic to a file;
otherwise the current device is used. Font and colour choices follow
[`theme_relab`](https://robustecologies.github.io/lucifer/reference/theme_relab.md).

## References

Durbin, J., & Koopman, S. J. (2012). *Time Series Analysis by State
Space Methods* (2nd ed.). Oxford University Press. ISBN 9780199641178.

## See also

[`SSM`](https://robustecologies.github.io/lucifer/reference/SSM.md),
[`as.demonoid.ssm_fit`](https://robustecologies.github.io/lucifer/reference/as.demonoid.ssm_fit.md),
[`print.ssm_fit`](https://robustecologies.github.io/lucifer/reference/print.ssm_fit.md),
[`summary.ssm_fit`](https://robustecologies.github.io/lucifer/reference/summary.ssm_fit.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving plot.ssm_fit
} # }
```
