# Mix two colors

Blend two colors together at a specified ratio.

## Usage

``` r
mix_colors(color1, color2, ratio = 0.5)
```

## Arguments

- color1:

  Hex color code for the first color.

- color2:

  Hex color code for the second color.

- ratio:

  Numeric value between 0 and 1 specifying the blend ratio. 0 = pure
  color2, 1 = pure color1, 0.5 = equal mix (default).

## Value

A hex color code representing the mixed color.

## References

Smith, A. R. (1978). Color gamut transform pairs. *ACM SIGGRAPH Computer
Graphics*, 12(3), 12-19.
[doi:10.1145/965139.807361](https://doi.org/10.1145/965139.807361)

## Examples

``` r
mix_colors("#FF0000", "#0000FF")
#> [1] "#800080"
mix_colors("#FF0000", "#0000FF", ratio = 0.75)
#> [1] "#BF0040"
```
