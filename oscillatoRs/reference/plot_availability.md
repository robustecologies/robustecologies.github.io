# Plot data availability timeline

Creates a timeline showing the temporal coverage of each climate index,
colored by category.

## Usage

``` r
plot_availability(data, title = NULL)
```

## Arguments

- data:

  A tibble containing climate index data in long format. Must have
  columns `date`, `index`, and optionally `category`.

- title:

  Character. Plot title.

## Value

A ggplot2 object.

## References

Wilke, C. O. (2019). *Fundamentals of data visualization*. O'Reilly
Media. ISBN 978-1492031086

## Examples

``` r
if (FALSE) { # \dontrun{
data(climate_monthly)
plot_availability(climate_monthly)
} # }
```
