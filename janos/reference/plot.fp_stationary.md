# Plot method for fp_stationary objects

Visualizes the Fokker-Planck stationary density and, optionally, the
drift potential. The default plot shows the numerical stationary density
overlaid with the exact Boltzmann distribution when available.

## Usage

``` r
# S3 method for class 'fp_stationary'
plot(x, type = c("density", "potential", "both"), ...)
```

## Arguments

- x:

  An `fp_stationary` object.

- type:

  Character string selecting the plot type: `"density"` (default),
  `"potential"`, or `"both"`.

- ...:

  Additional arguments (ignored).

## Value

A ggplot2 object.
