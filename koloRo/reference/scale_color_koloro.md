# Discrete color scale using koloRo palettes

Apply a koloRo palette as a discrete color scale in ggplot2.

## Usage

``` r
scale_color_koloro(palette = "okabe_ito", direction = 1, ...)

scale_colour_koloro(palette = "okabe_ito", direction = 1, ...)
```

## Arguments

- palette:

  Character string specifying the palette name. Default is "okabe_ito".

- direction:

  Integer. If 1 (default), colors are used as-is. If -1, the order of
  colors is reversed.

- ...:

  Additional arguments passed to
  [`ggplot2::discrete_scale()`](https://ggplot2.tidyverse.org/reference/discrete_scale.html).

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

# Basic usage
ggplot(mtcars, aes(x = wt, y = mpg, color = factor(cyl))) +
  geom_point(size = 3) +
  scale_color_koloro()

# Use a specific palette
ggplot(mtcars, aes(x = wt, y = mpg, color = factor(cyl))) +
  geom_point(size = 3) +
  scale_color_koloro("tol_bright")

# Reverse color order
ggplot(mtcars, aes(x = wt, y = mpg, color = factor(cyl))) +
  geom_point(size = 3) +
  scale_color_koloro("viridis", direction = -1)
} # }
```
