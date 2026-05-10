# Summary method for chameleon_analysis objects

Generate a comprehensive summary of a chameleon attractor analysis,
including diagnostics from all pipeline stages: embedding, MEMD
decomposition, EVT metrics, and scale-dependent analysis.

## Usage

``` r
# S3 method for class 'chameleon_analysis'
summary(object, ...)
```

## Arguments

- object:

  Object of class `chameleon_analysis`.

- ...:

  Additional arguments (ignored).

## Value

An object of class `summary.chameleon_analysis` containing:

- params:

  List of analysis parameters.

- embedding_summary:

  Summary of the Takens embedding.

- memd_summary:

  Summary of the MEMD decomposition.

- evt_summary:

  Summary of the full-attractor EVT metrics.

- scale_summary:

  Summary of the scale-dependent metrics.

- is_chameleon:

  Detection results.

- statistical_test:

  Statistical test results (if performed).

## Details

This summary provides a multi-section report covering all stages of the
chameleon analysis pipeline. Each component summary can be accessed
individually for more detailed diagnostics.

The print method displays a formatted report with key findings from each
stage, along with the overall chameleon assessment and interpretation
guidance.

## See also

[`chameleon_analysis`](https://robustecologies.github.io/chamaeleon/reference/chameleon_analysis.md)
for running the analysis;
[`print.chameleon_analysis`](https://robustecologies.github.io/chamaeleon/reference/print.chameleon_analysis.md),
[`plot.chameleon_analysis`](https://robustecologies.github.io/chamaeleon/reference/plot.chameleon_analysis.md);
[`chameleon_test`](https://robustecologies.github.io/chamaeleon/reference/chameleon_test.md)
for statistical testing.
