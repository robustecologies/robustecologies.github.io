# Cumulative sample function

The Cumulative Sample Function (CSF) is a visual MCMC diagnostic in
which the user may select a measure (such as a variable, summary
statistic, or other diagnostic), and observe a plot of how the measure
changes over cumulative posterior samples from MCMC, such as the output
of
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).
This may be considered to be a generalized extension of the `cumuplot`
in the coda package, which is a more restrictive form of the cusum
diagnostic introduced by Yu and Myckland (1998).

Yu and Myckland (1998) suggest that CSF plots should be examined after
traditional trace plots seem convergent, and assert that faster mixing
chains (which are more desirable) result in CSF plots that are more
"hairy" (as opposed to smooth), though this is subjective and has been
debated.

## Usage

``` r
CSF(
  x,
  name,
  method = "Quantiles",
  quantiles = c(0.025, 0.5, 0.975),
  output = FALSE
)
```

## Arguments

- x:

  This is a vector of posterior samples from MCMC.

- name:

  This is an optional name for vector `x`, and is input as a quoted
  string, such as `name="theta"`.

- method:

  This is a measure that will be observed over the course of cumulative
  samples of `x`. It defaults to `method="Quantiles"`, and optional
  methods include: `"ESS"`, `"Geweke.Diagnostic"`, `"HPD"`,
  `"is.stationary"`, `"Kurtosis"`, `"MCSE"`, `"MCSE.bm"`, `"MCSE.sv"`,
  `"Mean"`, `"Mode"`, `"N.Modes"`, `"Precision"`, `"Quantiles"`, and
  `"Skewness"`.

- quantiles:

  This optional argument applies only when `method="Quantiles"`, in
  which case this vector indicates the probabilities that will be
  observed. It defaults to the median and 95% probability interval
  bounds (see
  [`p.interval`](https://robustecologies.github.io/lucifer/reference/p.interval.md)
  for more information).

- output:

  Logical. If `output=TRUE`, then the results of the measure over the
  course of the cumulative samples will be output as an object, either a
  vector or matrix, depending on the `method` argument. The `output`
  argument defaults to `FALSE`.

## Value

See Details.

## Details

When `method="ESS"`, the effective sample size (ESS) is observed as a
function of the cumulative samples of `x`. For more information, see the
[`ESS`](https://robustecologies.github.io/lucifer/reference/ESS.md)
function.

When `method="Geweke.Diagnostic"`, the Z-score output of the Geweke
diagnostic is observed as a function of the cumulative samples of `x`.
For more information, see the
[`Geweke.Diagnostic`](https://robustecologies.github.io/lucifer/reference/Geweke.Diagnostic.md)
function.

When `method="HPD"`, the Highest Posterior Density (HPD) interval is
observed as a function of the cumulative samples of `x`. For more
information, see the
[`p.interval`](https://robustecologies.github.io/lucifer/reference/p.interval.md)
function.

When `method="is.stationary"`, stationarity is logically tested and the
result is observed as a function of the cumulative samples of `x`. For
more information, see the
[`is.stationary`](https://robustecologies.github.io/lucifer/reference/is.stationary.md)
function.

When `method="Kurtosis"`, kurtosis is observed as a function of the
cumulative samples of `x`.

When `method="MCSE"`, the Monte Carlo Standard Error (MCSE) estimated
with the `IMPS` method is observed as a function of the cumulative
samples of `x`. For more information, see the
[`MCSE`](https://robustecologies.github.io/lucifer/reference/MCSE.md)
function.

When `method="MCSE.bm"`, the Monte Carlo Standard Error (MCSE) estimated
with the `batch.means` method is observed as a function of the
cumulative samples of `x`. For more information, see the
[`MCSE`](https://robustecologies.github.io/lucifer/reference/MCSE.md)
function.

When `method="MCSE.sv"`, the Monte Carlo Standard Error (MCSE) estimated
with the `sample.variance` method is observed as a function of the
cumulative samples of `x`. For more information, see the
[`MCSE`](https://robustecologies.github.io/lucifer/reference/MCSE.md)
function.

When `method="Mean"`, the mean is observed as a function of the
cumulative samples of `x`.

When `method="Mode"`, the estimated mode is observed as a function of
the cumulative samples of `x`. For more information, see the
[`Mode`](https://robustecologies.github.io/lucifer/reference/Mode.md)
function.

When `method="N.Modes"`, the estimated number of modes is observed as a
function of the cumulative samples of `x`. For more information, see the
[`Modes`](https://robustecologies.github.io/lucifer/reference/Mode.md)
function.

When `method="Precision"`, the precision (inverse variance) is observed
as a function of the cumulative samples of `x`.

When `method="Quantiles"`, the quantiles selected with the `quantiles`
argument are observed as a function of the cumulative samples of `x`.

When `method="Skewness"`, skewness is observed as a function of the
cumulative samples of `x`.

## References

Yu, B. and Myckland, P. (1997). "Looking at Markov Samplers through
Cusum Path Plots: A Simple Diagnostic Idea". *Statistics and Computing*,
8(3), p. 275–286.

## See also

[`burnin`](https://robustecologies.github.io/lucifer/reference/burnin.md),
[`ESS`](https://robustecologies.github.io/lucifer/reference/ESS.md),
[`Geweke.Diagnostic`](https://robustecologies.github.io/lucifer/reference/Geweke.Diagnostic.md),
[`is.stationary`](https://robustecologies.github.io/lucifer/reference/is.stationary.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`MCSE`](https://robustecologies.github.io/lucifer/reference/MCSE.md),
[`Mode`](https://robustecologies.github.io/lucifer/reference/Mode.md),
[`Modes`](https://robustecologies.github.io/lucifer/reference/Mode.md),
and
[`p.interval`](https://robustecologies.github.io/lucifer/reference/p.interval.md).

## Examples

``` r
if (FALSE) { # \dontrun{
library(lucifer)
x <- rnorm(1000)
CSF(x, method="ESS")
CSF(x, method="Geweke.Diagnostic")
CSF(x, method="HPD")
CSF(x, method="is.stationary")
CSF(x, method="Kurtosis")
CSF(x, method="MCSE")
CSF(x, method="MCSE.bm")
CSF(x, method="MCSE.sv")
CSF(x, method="Mean")
CSF(x, method="Mode")
CSF(x, method="N.Modes")
CSF(x, method="Precision")
CSF(x, method="Quantiles")
CSF(x, method="Skewness")
} # }
```
