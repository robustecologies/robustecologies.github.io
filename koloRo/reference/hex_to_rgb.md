# Convert hex color to RGB values

Converts hexadecimal color codes to RGB values (0-255 scale).

## Usage

``` r
hex_to_rgb(hex)
```

## Arguments

- hex:

  Character. Hexadecimal color code (with or without '#').

## Value

Numeric vector of length 3 containing R, G, B values (0-255).

## Examples

``` r
hex_to_rgb("#FF0000")
#> [1] 255   0   0
hex_to_rgb("00FF00")
#> [1]   0 255   0
```
