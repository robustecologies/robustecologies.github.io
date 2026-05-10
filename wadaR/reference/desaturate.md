# Desaturate a color

Reduce the saturation of a color, making it more gray.

## Usage

``` r
desaturate(color, amount = 0.5)
```

## Arguments

- color:

  Hex color code or vector of colors.

- amount:

  Amount to desaturate (0-1). 0 = no change, 1 = fully gray.

## Value

A hex color code or vector of hex color codes.

## References

Smith, A. R. (1978). Color gamut transform pairs. *ACM SIGGRAPH Computer
Graphics*, 12(3), 12-19.
[doi:10.1145/965139.807361](https://doi.org/10.1145/965139.807361)

## Examples

``` r
desaturate("#FF0000", 0.5)
#>   #FF0000 
#> "#FF8080" 
desaturate(c("#FF0000", "#00FF00", "#0000FF"), 0.3)
#>   #FF0000   #00FF00   #0000FF 
#> "#FF4D4D" "#4DFF4D" "#4D4DFF" 
```
