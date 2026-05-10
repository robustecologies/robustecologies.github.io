# Six-level taxonomic-focus palette

Returns a named character vector of hex codes keyed by the six values of
the `taxonomic_focus` factor (`mammals`, `birds`, `fishes`,
`amphibians`, `reptiles`, `mixed`). The palette is colour-blind safe and
consistent with the Dark2 family used across other ecological
visualisations in the literature.

## Usage

``` r
vt_palette()
```

## Value

Named character vector of length 6 with hex codes.

## Details

Use as `scale_colour_manual(values = vt_palette())` or
`scale_fill_manual(values = vt_palette())`. The level order matches the
phylogenetic order used in the manuscript (Fishes -\> Amphibians -\>
Reptiles -\> Birds -\> Mammals -\> Mixed) so that legends are stable
across figures.

## References

Brewer, C. A. (2003). *A transition in improving maps: the ColorBrewer
example*. Cartography and Geographic Information Science, 30(2),
159-162.
[doi:10.1559/152304003100011126](https://doi.org/10.1559/152304003100011126)
.

## See also

[`vt_theme()`](https://robustecologies.github.io/VerteTIME/reference/vt_theme.md),
[`vt_plot_world_map()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_world_map.md),
[`vt_plot_alpha_timeline()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_alpha_timeline.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(ggplot2)
df <- data.frame(
  taxonomic_focus = factor(c("fishes","amphibians","reptiles",
                             "birds","mammals","mixed"))
)
ggplot(df, aes(taxonomic_focus, fill = taxonomic_focus)) +
  geom_bar() +
  scale_fill_manual(values = vt_palette()) +
  vt_theme()
} # }
```
