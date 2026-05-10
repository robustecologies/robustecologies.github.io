# Summary method for takens_embedding objects

Compute detailed diagnostics for a Takens embedding, including
coordinate statistics, delay correlations, and filling ratio. These
diagnostics help assess the quality of the embedding parameters.

## Usage

``` r
# S3 method for class 'takens_embedding'
summary(object, ...)
```

## Arguments

- object:

  Object of class `takens_embedding`.

- ...:

  Additional arguments (ignored).

## Value

An object of class `summary.takens_embedding` containing:

- m:

  Embedding dimension.

- tau:

  Time delay.

- n_original:

  Original time series length.

- n_embedded:

  Number of embedded points.

- window_size:

  Embedding window in samples.

- coord_stats:

  Data frame with mean, sd, min, max for each coordinate.

- correlations:

  Matrix of pairwise correlations between delay coordinates.

- lag1_correlations:

  Vector of correlations between adjacent coordinates.

- filling_ratio:

  Ratio of embedded to original length, indicating data utilization.

## Details

The delay coordinate correlations provide insight into the
appropriateness of the time delay tau. Ideally, adjacent delay
coordinates should have moderate correlation (0.3-0.7); very high
correlations indicate tau is too small (redundant information), while
very low correlations may indicate tau is too large (loss of dynamical
connection).

The filling ratio indicates what fraction of the original time series
contributes to the embedded trajectory. Low filling ratios occur with
large m or tau values relative to series length.

## See also

[`takens_embed`](https://robustecologies.github.io/chamaeleon/reference/takens_embed.md)
for creating embeddings,
[`estimate_embedding_params`](https://robustecologies.github.io/chamaeleon/reference/estimate_embedding_params.md)
for parameter estimation.
