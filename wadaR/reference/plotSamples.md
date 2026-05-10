# Plot samples

This function provides basic plots that are extended to include samples.

## Usage

``` r
plotSamples(X, Style = "KDE", LB = 0.025, UB = 0.975, Title = NULL, col = NULL)
```

## Arguments

- X:

  A \\N \times S\\ numerical matrix of \\N\\ records and \\S\\ samples.

- Style:

  A quoted string accepting the following values: "barplot", "dotchart",
  "hist", "KDE", or "Time-Series". Defaults to `Style="KDE"`.

- LB:

  The lower bound of a probability interval, which must be in the
  interval \[0,0.5).

- UB:

  The upper bound of a probability interval, which must be in the
  interval (0.5,1\].

- Title:

  Defaults to `NULL`, and otherwise accepts a quoted string that will be
  the title of the plot.

- col:

  Optional character vector of colors. When `NULL` (default), the RElab
  palette is used.

## Value

See Details.

## Details

The `plotSamples` function extends several basic plots from points to
samples. For example, it is common to use the `hist` function to plot a
histogram from a column vector. However, the user may desire to plot a
histogram of a column vector that was sampled numerous times, rather
than a simple column vector, in which a (usually 95%) probability
interval is also plotted to show the uncertainty around the sampled
median of each bin in the histogram.

The `plotSamples` function extends the `barplot`, `dotchart`, and `hist`
functions to include uncertainty due to samples. The `KDE` style of plot
is added so that a probability interval is shown around a sampled kernel
density estimate of a distribution, and the `Time-Series` style of plot
is added so that a probability interval is shown around a sampled
univariate time-series.

For each style of plot, three quantiles are plotted: the lower bound
(LB), median, and upper bound (UB).

One of many potential Bayesian applications is to examine the
uncertainty in a predictive distribution.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
N <- 100
S <- 100
X <- matrix(rnorm(N*S),N,S)
rownames(X) <- 1:100
plotSamples(X, Style="barplot", LB=0.025, UB=0.975)
plotSamples(X[1:10,], Style="dotchart", LB=0.025, UB=0.975)
plotSamples(X, Style="hist", LB=0.025, UB=0.975)
plotSamples(X, Style="KDE", LB=0.025, UB=0.975)
plotSamples(X, Style="Time-Series", LB=0.025, UB=0.975)
} # }
```
