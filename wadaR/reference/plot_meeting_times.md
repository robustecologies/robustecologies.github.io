# Plot meeting-time distribution from multiple coupled MCMC runs

Produces a histogram of meeting times from a list of `coupled_mcmc`
objects. The meeting-time distribution characterizes the mixing rate of
the underlying Markov kernel: heavier tails indicate slower convergence.

## Usage

``` r
plot_meeting_times(fits, col = NULL)
```

## Arguments

- fits:

  A list of objects of class `coupled_mcmc`, typically produced by
  `replicate(n, coupled_mcmc(...), simplify = FALSE)`.

- col:

  Optional character vector of hex color strings. When non-`NULL`,
  overrides the default RElab contrasting palette.

## Value

Invisibly returns the ggplot object.

## Details

Implementation of `plot_meeting_times`. Refer to the package vignettes
and the cited references for a complete algorithmic and mathematical
description.

## See also

[`coupled_mcmc`](https://robustecologies.github.io/lucifer/reference/coupled_mcmc.md),
[`plot.coupled_mcmc`](https://robustecologies.github.io/lucifer/reference/plot.coupled_mcmc.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fits <- replicate(50, coupled_mcmc(Model, Data,
  Initial.Values = c(0, 0), max_iterations = 5000),
  simplify = FALSE)
plot_meeting_times(fits)
} # }
```
