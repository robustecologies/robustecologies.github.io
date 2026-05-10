# Plot method for Wada analysis results

Creates a multi-panel summary visualization showing results from all
Wada detection methods applied to the basin structure.

## Usage

``` r
# S3 method for class 'wada_analysis'
plot(x, caption = TRUE, ...)
```

## Arguments

- x:

  A wada_analysis object.

- caption:

  Logical. If TRUE (default), render a one-line caption with the
  function name on the outer patchwork annotation. Child panels never
  carry their own caption to avoid duplication under composition.

- ...:

  Additional arguments (ignored).

## Value

A ggplot2 patchwork object, invisibly.

## References

Daza, A., Wagemakers, A., Georgeot, B., Guery-Odelin, D., & Sanjuan, M.
A. F. (2016). Basin entropy: a new tool to analyze uncertainty in
dynamical systems. *Scientific Reports*, 6, 31416.
[doi:10.1038/srep31416](https://doi.org/10.1038/srep31416)
