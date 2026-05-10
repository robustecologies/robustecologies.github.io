# Plot multiple palettes in a grid

Display multiple palettes in a grid layout for quick visual comparison.
Each palette appears as a compact strip with its name as the title,
allowing efficient overview of many palettes simultaneously.

## Usage

``` r
plot_palette_grid(palette_names, ncol = 2)
```

## Arguments

- palette_names:

  Character vector of palette names.

- ncol:

  Number of columns in the grid layout (default 2).

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
# Display colorblind-safe palettes in a 2-column grid
plot_palette_grid(c("okabe_ito", "wong", "tol_bright", "tol_muted"))

# Display Alhambra palettes in a 3-column grid
plot_palette_grid(c(
  "alhambra_nazari", "alhambra_zellige", "alhambra_blues",
  "alhambra_greens", "alhambra_reds", "alhambra_golds"
), ncol = 3)
} # }
```
