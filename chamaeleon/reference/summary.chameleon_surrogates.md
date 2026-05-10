# Summary method for chameleon_surrogates objects

Marginal statistics across a surrogate ensemble produced by
[`generate_surrogates`](https://robustecologies.github.io/chamaeleon/reference/generate_surrogates.md).
For `input_type == "scale_metrics"` the per-scale mean and standard
deviation of the dimension and persistence metrics are returned; for
`input_type == "series"` the per-element mean and standard deviation
across surrogate realisations are returned.

## Usage

``` r
# S3 method for class 'chameleon_surrogates'
summary(object, ...)
```

## Arguments

- object:

  Object of class `chameleon_surrogates`.

- ...:

  Additional arguments (ignored).

## Value

An object of class `summary.chameleon_surrogates` with components
`method`, `n_surrogates`, `input_type`, and either `D_mean`, `D_sd`,
`theta_mean`, `theta_sd`, `frequencies` (scale_metrics input) or
`value_mean`, `value_sd` (series input).

## See also

[`generate_surrogates`](https://robustecologies.github.io/chamaeleon/reference/generate_surrogates.md),
[`print.chameleon_surrogates`](https://robustecologies.github.io/chamaeleon/reference/print.chameleon_surrogates.md),
[`plot.chameleon_surrogates`](https://robustecologies.github.io/chamaeleon/reference/generate_surrogates.md).
