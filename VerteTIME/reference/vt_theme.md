# VerteTIME ggplot2 theme

Project-wide ggplot2 theme used by every figure shipped with the package
and by every plot method of the package S3 classes. Stacks
[`ggplot2::theme_minimal`](https://ggplot2.tidyverse.org/reference/ggtheme.html)
with bold title, grey-25 subtitle, grey-40 left-aligned size-8 caption
(per the project plotting convention) and a clean grid.

## Usage

``` r
vt_theme(base_size = 11)
```

## Arguments

- base_size:

  Base font size in points. Defaults to 11.

## Value

A `theme` object suitable for the `+` operator with any ggplot2 plot.
Re-applying overrides previous theme settings.

## Details

Every package plot calls `+ vt_theme()` after the geom layer is set up.
The subtitle slot carries one short technical sentence naming the
quantity shown; the caption slot is grey, left-aligned, and explains the
data subset, source and any non-default convention. This is enforced
project-wide so the rendered manuscript and the rendered website share a
consistent visual register.

## References

Wickham, H. (2016). *ggplot2: Elegant Graphics for Data Analysis*.
Springer-Verlag.
[doi:10.1007/978-3-319-24277-4](https://doi.org/10.1007/978-3-319-24277-4)
.

## See also

[`vt_palette()`](https://robustecologies.github.io/VerteTIME/reference/vt_palette.md),
[`vt_plot_whittaker()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_whittaker.md),
[`vt_plot_world_map()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_world_map.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(ggplot2)
ggplot(mtcars, aes(wt, mpg)) +
  geom_point() +
  labs(title = "Demo", subtitle = "fuel economy vs weight",
       caption = "mtcars dataset; for theme-only demonstration") +
  vt_theme()
} # }
```
