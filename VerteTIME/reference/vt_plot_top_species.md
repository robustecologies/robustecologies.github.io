# Top-N most monitored species across the compilation

Horizontal bar chart of the `n` species with the largest number of
observation rows in the compilation, coloured by the number of distinct
datasets that record them. Reads `species_id` from `x$observations`.

## Usage

``` r
vt_plot_top_species(x, n = 20)
```

## Arguments

- x:

  A `vt_compilation`.

- n:

  Integer; default 20. Number of species to list.

## Value

A `ggplot` object: horizontal bar chart with the species ordered from
most to least observed.

## References

Brlík, V., Šilarová, E., Škorpilová, J., Alonso, H., Anton, M., Aunins,
A., Benkö, Z., Biver, G., Busch, M., Chodkiewicz, T., et al. (2021).
*Long-term and large-scale multispecies dataset tracking population
changes of common European breeding birds*. Scientific Data, 8, 21.
[doi:10.1038/s41597-021-00804-2](https://doi.org/10.1038/s41597-021-00804-2)
.

Dornelas, M., Gotelli, N. J., McGill, B., Shimadzu, H., Moyes, F.,
Sievers, C., & Magurran, A. E. (2014). *Assemblage time series reveal
biodiversity change but not systematic loss*. Science, 344(6181),
296-299.
[doi:10.1126/science.1248484](https://doi.org/10.1126/science.1248484) .

## See also

[`vt_plot_top_genera()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_top_genera.md),
[`vt_plot_focus_bars()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_focus_bars.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(vertetime)
vt_plot_top_species(vertetime, n = 25)
} # }
```
