# Diverging color scale using koloRo palettes

Apply a koloRo diverging palette as a color scale in ggplot2. Ideal for
data with a meaningful midpoint (e.g., zero, average).

## Usage

``` r
scale_color_koloro_div(palette = "vik", midpoint = 0, direction = 1, ...)

scale_colour_koloro_div(palette = "vik", midpoint = 0, direction = 1, ...)
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

# Data with positive and negative values
df <- data.frame(
  x = 1:100,
  y = rnorm(100)
)

ggplot(df, aes(x, y, color = y)) +
  geom_point(size = 3) +
  scale_color_koloro_div()

# Use brown-blue-green diverging palette
ggplot(df, aes(x, y, color = y)) +
  geom_point(size = 3) +
  scale_color_koloro_div("brown_blue_green")
} # }
```
