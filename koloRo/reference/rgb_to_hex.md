# Convert RGB values to hex color

Converts RGB values (0-255 scale) to hexadecimal color code.

## Usage

``` r
rgb_to_hex(r, g, b)
```

## Arguments

- r:

  Numeric. Red value (0-255).

- g:

  Numeric. Green value (0-255).

- b:

  Numeric. Blue value (0-255).

## Value

Character. Hexadecimal color code with '#' prefix.

## Examples

``` r
rgb_to_hex(255, 0, 0)
#> [1] "#FF0000"
rgb_to_hex(0, 128, 255)
#> [1] "#0080FF"
```
