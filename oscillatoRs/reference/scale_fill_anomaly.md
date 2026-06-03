# Anomaly fill scale

Color scale for climate anomalies with blue for negative values, white
at zero, and red for positive values.

## Usage

``` r
scale_fill_anomaly(...)
```

## Arguments

- ...:

  Additional arguments passed to scale_fill_gradient2.

## Value

A ggplot2 scale object.

## References

Wilke, C. O. (2019). *Fundamentals of data visualization*. O'Reilly
Media. ISBN 978-1492031086

## Examples

``` r
if (FALSE) { # \dontrun{
library(ggplot2)
data(climate_monthly)
nao <- climate_monthly[climate_monthly$index == "NAO", ]

ggplot(nao, aes(x = date, y = value)) +
  geom_col(aes(fill = value)) +
  scale_fill_anomaly() +
  theme_climate()
} # }
```
