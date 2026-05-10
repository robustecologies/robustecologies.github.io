# Continuous color scale using koloRo palettes

Apply a koloRo palette as a continuous color gradient in ggplot2. Useful
for representing continuous numeric data.

## Usage

``` r
scale_color_koloro_c(palette = "viridis", direction = 1, ...)

scale_colour_koloro_c(palette = "viridis", direction = 1, ...)
```

## Arguments

- palette:

  Character string specifying the palette name. Default is "viridis".

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

# Basic usage with viridis
ggplot(mtcars, aes(x = wt, y = mpg, color = hp)) +
  geom_point(size = 3) +
  scale_color_koloro_c()

# Use plasma palette
ggplot(mtcars, aes(x = wt, y = mpg, color = hp)) +
  geom_point(size = 3) +
  scale_color_koloro_c("plasma")

# Use Alhambra sunset palette
ggplot(mtcars, aes(x = wt, y = mpg, color = hp)) +
  geom_point(size = 3) +
  scale_color_koloro_c("alhambra_sunset")
} # }
```
