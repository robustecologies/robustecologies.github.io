# Climate visualization theme

A clean, publication-ready ggplot2 theme for climate index
visualizations.

## Usage

``` r
theme_climate(base_size = 11, base_family = "")
```

## Arguments

- base_size:

  Numeric. Base font size in points.

- base_family:

  Character. Base font family.

## Value

A ggplot2 theme object.

## References

Wilke, C. O. (2019). *Fundamentals of data visualization*. O'Reilly
Media. ISBN 978-1492031086

## Examples

``` r
library(ggplot2)
ggplot(mtcars, aes(mpg, wt)) +
  geom_point() +
  theme_climate()
```
