# Plot samples from the output of Lucifer

This may be used to plot, or save plots of, samples in an object of
class `demonoid`. The `type` argument selects among several
visualization styles including trace plots, density plots,
autocorrelation, interval summaries, algorithm-specific diagnostics, and
more.

## Usage

``` r
# S3 method for class 'demonoid'
plot(
  x,
  BurnIn = 0,
  Data = NULL,
  PDF = FALSE,
  Parms = NULL,
  FileName = paste0("laplacesDemon-plot_", format(Sys.time(), "%Y-%m-%d_%T"), ".pdf"),
  type = "diagnostics",
  ground_truth = NULL,
  col = NULL,
  ...
)
```

## Arguments

- x:

  This required argument is an object of class `demonoid`.

- BurnIn:

  This argument requires zero or a positive integer that indicates the
  number of thinned samples to discard as burn-in for the purposes of
  plotting.

- Data:

  Optional. The list of data supplied to
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).
  Monitor variable names are recovered automatically from the fit
  object, so `Data` is no longer required for standard plot types.

- PDF:

  This logical argument indicates whether or not the user wants Lucifer
  to save the plots as a .pdf file.

- Parms:

  This argument accepts a vector of quoted strings to be matched for
  selecting parameters for plotting. This argument defaults to `NULL`
  and selects every parameter for plotting.

- FileName:

  This argument accepts a string and saves the plot under the specified
  name.

- type:

  Character. Plot type to produce. One of `"diagnostics"` (default),
  `"trace"`, `"density"`, `"histogram"`, `"acf"`, `"all"`,
  `"intervals"`, `"areas"`, `"pairs"`, `"rank"`, `"rhat"`, `"neff"`,
  `"energy"`, `"divergences"`, `"nuts"`, or `"parcoord"`. The `"all"`
  type shows trace, density, and ACF in a 3-column grid with one
  parameter per row.

- ground_truth:

  Optional named numeric vector of true parameter values. When provided,
  reference lines or markers are overlaid on the appropriate plot
  elements.

- col:

  Optional character vector of hex color strings. When non-`NULL`,
  overrides the default RElab contrasting palette for all plot elements.
  The first color becomes the primary, the second the reference, and so
  on; for multi-chain plots the colors are used in order.

- ...:

  Additional arguments passed to sub-plot functions.

## Value

Invisibly returns `x`. Called for its graphical side effect.

## Details

When `type = "diagnostics"` (the default), plots are arranged in a \\3
\times 3\\ matrix. Each row represents a parameter, the deviance, or a
monitored variable. The left column displays trace plots, the middle
column displays kernel density plots, and the right column displays
autocorrelation (ACF) plots.

Additional types dispatch to specialized visualization functions:
`"trace"` and `"density"` show faceted single-panel views with optional
`ground_truth` overlays; `"acf"` shows faceted ACF plots; `"intervals"`
delegates to
[`caterpillar.plot`](https://robustecologies.github.io/lucifer/reference/caterpillar.plot.md);
`"areas"` shows ridgeline density areas; `"pairs"`, `"rank"`, `"rhat"`,
`"neff"`, `"energy"`, `"divergences"`, `"nuts"`, and `"parcoord"`
dispatch to the algorithm-aware MCMC diagnostic functions in
`plot_mcmc.R`.

## See also

[`burnin`](https://robustecologies.github.io/lucifer/reference/burnin.md),
[`ESS`](https://robustecologies.github.io/lucifer/reference/ESS.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`caterpillar.plot`](https://robustecologies.github.io/lucifer/reference/caterpillar.plot.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Define model and data
N <- 25
y <- rnorm(N, 2, 0.5)
Data <- list(N = N, y = y, mon.names = "LP",
  parm.names = c("mu", "log.sigma"))
Model <- function(parm, Data) {
  mu <- parm[1]
  sigma <- exp(parm[2])
  LL <- sum(dnorm(Data$y, mu, sigma, log = TRUE))
  LP <- LL + dnormv(mu, 0, 1000, log = TRUE) +
    dhalfcauchy(sigma, 25, log = TRUE)
  yhat <- rnorm(Data$N, mu, sigma)
  Monitor <- LP
  return(list(LP = LP, Dev = -2 * LL, Monitor = Monitor,
    yhat = yhat, parm = parm))
}

# Fit model
fit <- lucifer(Model, Data, Initial.Values = c(0, 0),
  Iterations = 1000, Status = 200, Thinning = 1,
  Algorithm = "CHARM", Specs = NULL)

# Plot all diagnostics (Data no longer required)
plot(fit)

# Trace plots with ground truth
plot(fit, type = "trace",
  ground_truth = c(mu = 2, log.sigma = log(0.5)))

# Density plots
plot(fit, type = "density")

# Caterpillar / interval plot
plot(fit, type = "intervals")

# Pairs plot
plot(fit, type = "pairs")
} # }
```
