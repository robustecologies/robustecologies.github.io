# Create delay embedding from time series

Constructs a delay embedding of a scalar time series using either
standard delay coordinates or peak embedding (von Kármán style).

## Usage

``` r
turbulence_embed(
  x,
  embed_dim = 3L,
  tau = 10L,
  method = c("delay", "peaks"),
  min_separation = 0.1,
  smooth_window = 0
)
```

## Arguments

- x:

  Time series (numeric vector) or simulation object.

- embed_dim:

  Embedding dimension (default 3).

- tau:

  Time delay for embedding (in samples).

- method:

  Embedding method: "delay" for standard Takens embedding, "peaks" for
  peak extraction (requires simulation object).

- min_separation:

  Minimum time separation for peak extraction (only used if
  method="peaks").

- smooth_window:

  Smoothing window for peak detection.

## Value

Object of class "turbulence_embedding" containing the embedded
coordinates and metadata.

## Details

Standard delay embedding reconstructs the attractor using coordinates:
\$\$(x(t), x(t-\tau), x(t-2\tau), \ldots, x(t-(m-1)\tau))\$\$

Peak embedding extracts local maxima with minimum time separation, then
creates embedding vectors from consecutive peaks. This is particularly
effective for oscillatory systems.

## References

Takens, F. (1981). Detecting strange attractors in turbulence. *Lecture
Notes in Mathematics*, 898, 366-381.
[doi:10.1007/BFb0091924](https://doi.org/10.1007/BFb0091924)

## See also

[`print.turbulence_embedding()`](https://robustecologies.github.io/tuRbulence/reference/print.turbulence_embedding.md),
[`summary.turbulence_embedding()`](https://robustecologies.github.io/tuRbulence/reference/summary.turbulence_embedding.md),
[`plot.turbulence_embedding()`](https://robustecologies.github.io/tuRbulence/reference/plot.turbulence_embedding.md),
[`turbulence_lyapunov()`](https://robustecologies.github.io/tuRbulence/reference/turbulence_lyapunov.md),
[`turbulence_cao()`](https://robustecologies.github.io/tuRbulence/reference/turbulence_cao.md),
[`turbulence_sim()`](https://robustecologies.github.io/tuRbulence/reference/turbulence_sim.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate and embed
sim <- vonkarman_sim(mu = 0.3, n_steps = 100000, seed = 42)

# Delay embedding
emb_delay <- turbulence_embed(sim, embed_dim = 3, tau = 10, method = "delay")
print(emb_delay)

# Peak embedding (von Karman style)
emb_peaks <- turbulence_embed(sim, embed_dim = 3, method = "peaks")
print(emb_peaks)
} # }
```
