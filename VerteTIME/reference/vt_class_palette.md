# Class-level colour palette aligned with the GBIF backbone classes

Returns a named character vector of hex codes keyed by the most common
GBIF classes that surface in the VerteTIME compilation (`Aves`,
`Mammalia`, `Squamata`, `Testudines`, `Amphibia`, `Actinopterygii`,
`Elasmobranchii`, `Insecta`, `Gastropoda`, `Other`). The Reptilia
families (`Squamata`, `Testudines`) share a hue; the fish classes share
a hue; non-vertebrate classes that surface in the compilation by way of
the species-name regex (e.g. `Insecta`, `Gastropoda`) get a single grey
to make them visually subordinate.

## Usage

``` r
vt_class_palette()
```

## Value

Named character vector of hex codes.

## Details

Use as
`scale_colour_manual(values = vt_class_palette(), na.value = "grey60")`.
The class-level palette is a thin wrapper of
`RColorBrewer::brewer.pal(8, "Dark2")` so figures using both
[`vt_palette()`](https://robustecologies.github.io/VerteTIME/reference/vt_palette.md)
and `vt_class_palette()` read consistently.

## References

Brewer, C. A. (2003). *A transition in improving maps: the ColorBrewer
example*. Cartography and Geographic Information Science, 30(2),
159-162.
[doi:10.1559/152304003100011126](https://doi.org/10.1559/152304003100011126)
.

## See also

[`vt_palette()`](https://robustecologies.github.io/VerteTIME/reference/vt_palette.md),
[`vt_theme()`](https://robustecologies.github.io/VerteTIME/reference/vt_theme.md)

## Examples

``` r
if (FALSE) { # \dontrun{
vt_class_palette()
} # }
```
