# Diverging fill scale using koloRo palettes

Apply a koloRo diverging palette as a fill scale in ggplot2.

## Usage

``` r
scale_fill_koloro_div(palette = "vik", midpoint = 0, direction = 1, ...)
```

## Arguments

- palette:

  Character string specifying the palette name. Default is "vik"
  (blue-white-red diverging).

- midpoint:

  Numeric value specifying the data value that should correspond to the
  middle color of the palette. Default is 0.

- direction:

  Integer. If 1 (default), colors are used as-is. If -1, the order of
  colors is reversed.

- ...:

  Additional arguments passed to
  [`ggplot2::scale_color_gradientn()`](https://ggplot2.tidyverse.org/reference/scale_gradient.html).

## Value

A ggplot2 scale object.

## References

Wickham, H. (2016). *ggplot2: Elegant graphics for data analysis* (2nd
ed.). Springer-Verlag.
[doi:10.1007/978-3-319-24277-4](https://doi.org/10.1007/978-3-319-24277-4)

## Examples

``` r
if (FALSE) { # \dontrun{
library(ggplot2)

# Heatmap with diverging colors
df <- expand.grid(x = 1:10, y = 1:10)
df$z <- rnorm(100)

ggplot(df, aes(x, y, fill = z)) +
  geom_tile() +
  scale_fill_koloro_div()
} # }
```
