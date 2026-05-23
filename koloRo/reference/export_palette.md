# Export palette to various formats

Export a palette to different color formats (hex, RGB, HSL) and
optionally save to a file.

## Usage

``` r
export_palette(
  palette_name,
  format = c("hex", "rgb", "rgb_normalized"),
  file = NULL
)
```

## Arguments

- palette_name:

  Character string specifying the palette name.

- format:

  Character string specifying the output format:

  - `"hex"`: Hexadecimal color codes (default)

  - `"rgb"`: RGB values (0-255)

  - `"rgb_normalized"`: RGB values (0-1)

- file:

  Optional character string specifying a file path to save the output.
  If NULL (default), returns the result without saving.

## Value

A data frame with the palette colors in the specified format.

## References

Itten, J. (1961). *The art of color*. Reinhold Publishing. ISBN
978-0471289296

CIE (2004). Colorimetry (3rd ed.). CIE Publication 15:2004. Commission
Internationale de l'Eclairage. ISBN 978-3901906336

## Examples

``` r
if (FALSE) { # \dontrun{
# Get hex colors
export_palette("viridis", format = "hex")

# Get RGB values
export_palette("okabe_ito", format = "rgb")

# Save to CSV
export_palette("tol_bright", format = "rgb", file = "tol_bright.csv")
} # }
```
