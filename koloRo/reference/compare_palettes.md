# Compare multiple palettes

Creates a side-by-side visual comparison of multiple palettes using base
R graphics. Each palette is displayed as a horizontal strip with the
palette name on the left, allowing direct comparison of color schemes.

## Usage

``` r
compare_palettes(palette_names, show_hex = FALSE)
```

## Arguments

- palette_names:

  Character vector of palette names to compare.

- show_hex:

  Logical. If TRUE, displays hex codes below colors when palettes have
  10 or fewer colors. Default is FALSE.

## Value

Invisibly returns NULL. Called for its side effect of producing a plot.

## Details

The function automatically adjusts the plot layout based on the number
of palettes being compared. For more flexible visualization options,
consider using
[`plot_palette`](https://robustecologies.github.io/koloRo/reference/plot_palette.md)
with a vector of palette names.

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
# Compare gold standard colorblind-safe palettes
compare_palettes(c("okabe_ito", "wong", "tol_bright"))

# Compare Alhambra location palettes
compare_palettes(c(
  "alhambra_patio_leones",
  "alhambra_sala_abencerrajes",
  "alhambra_generalife"
))

# Compare with hex codes shown
compare_palettes(c("viridis", "magma", "plasma"), show_hex = TRUE)
} # }
```
