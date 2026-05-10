# List available color palettes

Returns a data frame summarizing all available color palettes or those
within a specific category. This function is useful for discovering
palette options and understanding the scope of the koloRo collection.

## Usage

``` r
list_palettes(category = NULL)
```

## Arguments

- category:

  Character string specifying the palette category to list. If NULL
  (default), returns all palettes from all categories.

## Value

A data frame with columns:

- `palette_name`: Name of the palette (use with
  [`palettes()`](https://robustecologies.github.io/koloRo/reference/palettes.md))

- `category`: Category the palette belongs to

- `n_colors`: Number of colors in the palette

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
# List all available palettes
all_palettes <- list_palettes()
head(all_palettes)

# List only colorblind-safe palettes
cb_palettes <- list_palettes(category = "colorblind")
print(cb_palettes)

# List Alhambra historical palettes
alhambra_list <- list_palettes(category = "alhambra")
nrow(alhambra_list)  # Count available palettes
} # }
```
