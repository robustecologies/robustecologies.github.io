# Binned color scale using koloRo palettes

Apply a koloRo palette as a binned color scale in ggplot2. Useful for
creating stepped color gradients.

## Usage

``` r
scale_color_koloro_binned(palette = "viridis", n_bins = 6, direction = 1, ...)

scale_colour_koloro_binned(palette = "viridis", n_bins = 6, direction = 1, ...)
```

## Arguments

- palette:

  Character string specifying the palette name. Default is "viridis".

- n_bins:

  Integer specifying the number of bins. Default is 6.

- direction:

  Integer. If 1 (default), colors are used as-is. If -1, the order of
  colors is reversed.

- ...:

  Additional arguments passed to
  [`ggplot2::scale_color_stepsn()`](https://ggplot2.tidyverse.org/reference/scale_steps.html).

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

ggplot(mtcars, aes(x = wt, y = mpg, color = hp)) +
  geom_point(size = 3) +
  scale_color_koloro_binned()

# More bins
ggplot(mtcars, aes(x = wt, y = mpg, color = hp)) +
  geom_point(size = 3) +
  scale_color_koloro_binned("plasma", n_bins = 10)
} # }
```
