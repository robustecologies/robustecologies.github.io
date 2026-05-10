# Plot all palettes from a category

Display all palettes belonging to a specific category in a stacked
comparison view. This provides a comprehensive overview of the available
options within a category.

## Usage

``` r
plot_category(category, max_display = 20)
```

## Arguments

- category:

  Character string specifying the category name (e.g., "alhambra",
  "colorblind", "scientific").

- max_display:

  Integer specifying the maximum number of palettes to display (default
  20). If the category contains more palettes, a message indicates how
  many are shown.

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
# View all colorblind-safe palettes
plot_category("colorblind")

# View Alhambra palettes (limited to first 15)
plot_category("alhambra", max_display = 15)

# View scientific visualization palettes
plot_category("scientific")
} # }
```
