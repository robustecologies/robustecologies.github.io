# Get complementary color

Calculate the complementary color (opposite on the color wheel) for a
given color.

## Usage

``` r
complementary_color(color)
```

## Arguments

- color:

  Hex color code.

## Value

A hex color code representing the complementary color.

## References

Smith, A. R. (1978). Color gamut transform pairs. *ACM SIGGRAPH Computer
Graphics*, 12(3), 12-19.
[doi:10.1145/965139.807361](https://doi.org/10.1145/965139.807361)

Itten, J. (1961). *The art of color*. Reinhold Publishing. ISBN
978-0471289296

## Examples

``` r
complementary_color("#FF0000")
#> [1] "#00FFFF"
complementary_color("#1E4D8C")
#> [1] "#8C5D1E"
```
