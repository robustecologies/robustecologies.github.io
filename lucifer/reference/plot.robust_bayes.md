# Plot robust Bayesian sensitivity analysis results

Produces diagnostic visualizations for `robust_bayes` objects. The
default `"dashboard"` type provides a multi-panel overview of all
available modules. Individual plot types focus on specific analyses.

## Usage

``` r
# S3 method for class 'robust_bayes'
plot(x, type = "dashboard", Parms = NULL, col = NULL, ...)
```

## Arguments

- x:

  An object of class `robust_bayes`.

- type:

  Character string specifying the plot type. One of `"dashboard"`
  (default, multi-panel overview), `"power"` (power-scaling heatmap),
  `"power_density"` (overlaid density curves at representative alpha
  values, similar to priorsense), `"divergence"` (pairwise divergence
  tile plot), `"influence"` (Pareto-k and mean shift scatter),
  `"conflict"` (prior vs posterior forest plot), `"weights"` (grouped
  bar chart of model weights), or `"posterior"` (ridgeline densities
  across models).

- Parms:

  Character string or regex for subsetting parameters. Applies to
  `"power"`, `"conflict"`, and `"posterior"` plot types.

- col:

  Optional color vector override.

- ...:

  Additional arguments (currently ignored).

## Value

Invisibly returns the ggplot object(s).

## Details

Produces diagnostic plots of a RobustBayes prior-sensitivity analysis
produced by
[`RobustBayes`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md).
Summary of the content is given below. Default output renders a
multi-panel graphic (trace, density, and autocorrelation where
applicable). The `PDF` argument captures the graphic to a file;
otherwise the current device is used. Font and colour choices follow
[`theme_relab`](https://robustecologies.github.io/lucifer/reference/theme_relab.md).

## References

Berger, J. O., & Berliner, L. M. (1986). Robust Bayes and empirical
Bayes analysis with epsilon-contaminated priors. *Annals of Statistics*,
14(2), 461-486.
[doi:10.1214/aos/1176349933](https://doi.org/10.1214/aos/1176349933)

## See also

[`RobustBayes`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md)

## Examples

``` r
if (FALSE) { # \dontrun{
rb <- RobustBayes(fit, Model = Model, Data = MyData)
plot(rb)
plot(rb, type = "power")
plot(rb, type = "power_density")
plot(rb, type = "power_density", n_alpha = 9)
plot(rb, type = "influence")
} # }
```
