# Get detailed information about a palette

Returns comprehensive information about a specific palette including its
category, number of colors, hex values, and RGB values.

## Usage

``` r
palette_info(palette_name)
```

## Arguments

- palette_name:

  Character string specifying the palette name.

## Value

A list containing:

- `name`: Palette name

- `category`: Category the palette belongs to

- `n_colors`: Number of colors in the palette

- `colors`: Character vector of hex color codes

- `rgb`: Matrix of RGB values (one row per color)

## References

Itten, J. (1961). *The art of color*. Reinhold Publishing. ISBN
978-0471289296

CIE (2004). Colorimetry (3rd ed.). CIE Publication 15:2004. Commission
Internationale de l'Eclairage. ISBN 978-3901906336

## Examples

``` r
if (FALSE) { # \dontrun{
info <- palette_info("okabe_ito")
info$colors
info$rgb
} # }
```
