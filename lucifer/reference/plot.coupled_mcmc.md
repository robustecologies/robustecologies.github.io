# Plot method for coupled MCMC objects

Produces diagnostic plots for a single coupled MCMC run. The default
type shows trace plots of both chains with the meeting time marked. For
plotting meeting-time distributions from multiple runs, use
[`plot_meeting_times`](https://robustecologies.github.io/lucifer/reference/plot_meeting_times.md).

## Usage

``` r
# S3 method for class 'coupled_mcmc'
plot(x, type = NULL, col = NULL, ...)
```

## Arguments

- x:

  An object of class `coupled_mcmc`.

- type:

  Character: `"trace"` (default), `"meeting"`, or `"log_posterior"`.

- col:

  Optional character vector of hex color strings. When non-`NULL`,
  overrides the default RElab contrasting palette.

- ...:

  Additional arguments (currently ignored).

## Value

Invisibly returns the ggplot object(s).

## Details

Implementation of `plot.coupled_mcmc`. Refer to the package vignettes
and the cited references for a complete algorithmic and mathematical
description.

## See also

[`plot_meeting_times`](https://robustecologies.github.io/lucifer/reference/plot_meeting_times.md)
for histograms from multiple runs.

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving plot.coupled_mcmc
} # }
```
