# Get analogous colors

Generate analogous colors (adjacent on the color wheel) for a given
color.

## Usage

``` r
analogous_colors(color, n = 3, angle = 30)
```

## Arguments

- color:

  Hex color code.

- n:

  Number of analogous colors to generate (default 3, including the
  original color).

- angle:

  Angle in degrees between adjacent colors (default 30).

## Value

A character vector of hex color codes.

## References

Smith, A. R. (1978). Color gamut transform pairs. *ACM SIGGRAPH Computer
Graphics*, 12(3), 12-19.
[doi:10.1145/965139.807361](https://doi.org/10.1145/965139.807361)

Itten, J. (1961). *The art of color*. Reinhold Publishing. ISBN
978-0471289296

## Examples

``` r
analogous_colors("#FF0000")
#> [1] "#FF0080" "#FF0000" "#FF8000"
analogous_colors("#1E4D8C", n = 5, angle = 15)
#> [1] "#1E848C" "#1E698C" "#1E4D8C" "#1E318C" "#261E8C"
```
