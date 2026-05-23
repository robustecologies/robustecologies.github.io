# Add transparency to a color

Adds an alpha (transparency) channel to one or more colors. This is
useful for creating semi-transparent layers in visualizations,
particularly when overlapping elements need to remain visible.

## Usage

``` r
add_alpha(color, alpha = 0.5)
```

## Arguments

- color:

  Hex color code or vector of hex color codes.

- alpha:

  Numeric value between 0 and 1 specifying the opacity level. 0 is fully
  transparent, 1 is fully opaque. Default is 0.5.

## Value

A character vector of hex color codes with alpha channel (8-character
hex format: \#RRGGBBAA), same length as the input `color` vector.

## Details

The function converts the input color to RGB, then creates a new color
specification including the alpha channel. The output format is
compatible with R graphics devices that support alpha transparency.

## References

Smith, A. R. (1978). Color gamut transform pairs. *ACM SIGGRAPH Computer
Graphics*, 12(3), 12-19.
[doi:10.1145/965139.807361](https://doi.org/10.1145/965139.807361)

## Examples

``` r
if (FALSE) { # \dontrun{
base <- "#1E4D8C"  # Alhambra blue

# Add 50% transparency
semi_transparent <- add_alpha(base, 0.5)

# Create an opacity gradient
opacity_scale <- add_alpha(base, seq(0.1, 1, by = 0.1))

# Apply to multiple colors
cols <- palettes(palette = "okabe_ito")
transparent_cols <- add_alpha(cols, 0.7)
} # }
```
