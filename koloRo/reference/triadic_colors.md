# Get triadic colors

Generate triadic colors (evenly spaced 120 degrees apart on the color
wheel).

## Usage

``` r
triadic_colors(color)
```

## Arguments

- color:

  Hex color code.

## Value

A character vector of 3 hex color codes.

## References

Smith, A. R. (1978). Color gamut transform pairs. *ACM SIGGRAPH Computer
Graphics*, 12(3), 12-19.
[doi:10.1145/965139.807361](https://doi.org/10.1145/965139.807361)

Itten, J. (1961). *The art of color*. Reinhold Publishing. ISBN
978-0471289296

## Examples

``` r
triadic_colors("#FF0000")
#> [1] "#FF0000" "#00FF00" "#0000FF"
triadic_colors("#1E4D8C")
#> [1] "#1E4D8C" "#8C1E4D" "#4D8C1E"
```
