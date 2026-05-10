# RElab ggplot2 theme

A clean, publication-quality ggplot2 theme based on `theme_minimal` with
restrained styling: bold title, subtle gridlines, bottom legend, and
comfortable margins.

## Usage

``` r
theme_relab(base_size = 11, base_family = "")
```

## Arguments

- base_size:

  Base font size in points. Default 11.

- base_family:

  Base font family. Default `""`.

## Value

A ggplot2 theme object.

## Details

Implementation of `theme_relab`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
library(ggplot2)
ggplot(mtcars, aes(wt, mpg)) + geom_point() + theme_relab()
} # }
```
