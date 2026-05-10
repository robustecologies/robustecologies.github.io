# Plot color palette or palette ramp

Visualize color palettes. Can display single palettes, multiple palettes
for comparison, or interpolated ramps.

## Usage

``` r
plot_palette(x, show_hex = NULL, title = NULL, cex_text = 0.7)
```

## Arguments

- x:

  Can be:

  - Character vector of hex colors

  - Character string with palette name

  - Character vector with multiple palette names (for comparison)

- show_hex:

  Logical. Show hex codes below colors (default TRUE for \<15 colors).

- title:

  Character string for plot title (auto-generated if NULL).

- cex_text:

  Numeric. Text size for hex codes (default 0.7).

## Value

Invisibly returns NULL. Called for its side effect of producing a plot.

## References

Itten, J. (1961). *The art of color*. Reinhold Publishing. ISBN
978-0471289296

CIE (2004). Colorimetry (3rd ed.). CIE Publication 15:2004. Commission
Internationale de l'Eclairage. ISBN 978-3901906336

## Examples

``` r
if (FALSE) { # \dontrun{
# Plot a single palette by name
plot_palette("alhambra_nazari")

# Plot a palette ramp
plot_palette(palette_ramp("okabe_ito", n = 50))

# Plot directly from palettes()
plot_palette(palettes(palette = "viridis"))

# Compare multiple palettes
plot_palette(c("alhambra_nazari", "okabe_ito", "tol_bright"))

# Without hex codes
plot_palette("alhambra_zellige", show_hex = FALSE)
} # }
```
