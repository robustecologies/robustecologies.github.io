# Display color palette visually

Creates a simple visual display of a color palette using base R
graphics. Each color is shown as a large square with optional hex code
labels below. This function provides a quick way to preview a palette
before using it in visualizations.

## Usage

``` r
display_palette(palette_name, show_hex = TRUE, show_category = TRUE)
```

## Arguments

- palette_name:

  Character string specifying the palette name.

- show_hex:

  Logical. If TRUE (default), displays hex codes below colors.

- show_category:

  Logical. If TRUE (default), shows the category in the title.

## Value

Invisibly returns NULL. Called for its side effect of producing a plot.

## Details

This function modifies graphical parameters temporarily and restores
them after plotting. For comparing multiple palettes, use
[`plot_palette`](https://robustecologies.github.io/koloRo/reference/plot_palette.md)
or
[`compare_palettes`](https://robustecologies.github.io/koloRo/reference/compare_palettes.md)
instead.

## References

Itten, J. (1961). *The art of color*. Reinhold Publishing. ISBN
978-0471289296

CIE (2004). Colorimetry (3rd ed.). CIE Publication 15:2004. Commission
Internationale de l'Eclairage. ISBN 978-3901906336

Smith, A. R. (1978). Color gamut transform pairs. *ACM SIGGRAPH Computer
Graphics*, 12(3), 12-19.
[doi:10.1145/965139.807361](https://doi.org/10.1145/965139.807361)

Okabe, M., & Ito, K. (2008). Color universal design (CUD): how to make
figures and presentations that are friendly to colorblind people.
*J\*Fly Data Depository for Drosophila Researchers*.
<https://jfly.uni-koeln.de/color/>

Wickham, H. (2016). *ggplot2: Elegant graphics for data analysis* (2nd
ed.). Springer-Verlag.
[doi:10.1007/978-3-319-24277-4](https://doi.org/10.1007/978-3-319-24277-4)

## Examples

``` r
if (FALSE) { # \dontrun{
# Display a single palette with hex codes
display_palette("okabe_ito")

# Display without hex codes for cleaner view
display_palette("viridis", show_hex = FALSE)

# Display without category information
display_palette("alhambra_nazari", show_category = FALSE)
} # }
```
