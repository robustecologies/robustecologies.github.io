# Invert a color

Invert a color by subtracting its RGB values from 255.

## Usage

``` r
invert_color(color)
```

## Arguments

- color:

  Hex color code or vector of colors.

## Value

A hex color code or vector of hex color codes.

## References

Smith, A. R. (1978). Color gamut transform pairs. *ACM SIGGRAPH Computer
Graphics*, 12(3), 12-19.
[doi:10.1145/965139.807361](https://doi.org/10.1145/965139.807361)

## Examples

``` r
invert_color("#FF0000")
#>   #FF0000 
#> "#00FFFF" 
invert_color(c("#FFFFFF", "#000000"))
#>   #FFFFFF   #000000 
#> "#000000" "#FFFFFF" 
```
