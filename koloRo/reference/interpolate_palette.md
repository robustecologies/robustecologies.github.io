# Interpolate palette to specified number of colors

Create a smooth interpolation of a palette to generate the specified
number of colors.

## Usage

``` r
interpolate_palette(palette_name, n, reverse = FALSE)
```

## Arguments

- palette_name:

  Character string specifying the palette name, or a character vector of
  hex colors.

- n:

  Integer specifying the number of colors to generate.

- reverse:

  Logical. If TRUE, reverses the palette order before interpolation
  (default FALSE).

## Value

A character vector of hex color codes.

## References

Itten, J. (1961). *The art of color*. Reinhold Publishing. ISBN
978-0471289296

CIE (2004). Colorimetry (3rd ed.). CIE Publication 15:2004. Commission
Internationale de l'Eclairage. ISBN 978-3901906336

## Examples

``` r
if (FALSE) { # \dontrun{
# Interpolate viridis to 256 colors
cols <- interpolate_palette("viridis", n = 256)

# Interpolate custom colors
cols <- interpolate_palette(c("#FF0000", "#00FF00", "#0000FF"), n = 50)
} # }
```
