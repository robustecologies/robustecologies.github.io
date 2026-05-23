# Plot method for turbulence embedding objects

Creates a 2D or 3D attractor visualization from the delay or peak
embedding, with points colored by sequence index.

## Usage

``` r
# S3 method for class 'turbulence_embedding'
plot(x, type = c("2d", "3d"), dims = NULL, ...)
```

## Arguments

- x:

  A turbulence_embedding object.

- type:

  Character. "2d" for ggplot2 projection or "3d" for plotly (only if
  embed_dim \>= 3). Default is "2d".

- dims:

  Integer vector of length 2 or 3 specifying which embedding dimensions
  to plot. Default is c(1, 2) for 2D.

- ...:

  Additional arguments (ignored).

## Value

A ggplot2 or plotly object.

## See also

[`turbulence_embed()`](https://robustecologies.github.io/tuRbulence/reference/turbulence_embed.md),
[`print.turbulence_embedding()`](https://robustecologies.github.io/tuRbulence/reference/print.turbulence_embedding.md),
[`summary.turbulence_embedding()`](https://robustecologies.github.io/tuRbulence/reference/summary.turbulence_embedding.md),
[`turbulence_lyapunov()`](https://robustecologies.github.io/tuRbulence/reference/turbulence_lyapunov.md)
