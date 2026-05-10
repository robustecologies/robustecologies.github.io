# Alpha-diversity timeline

Trajectory of an alpha-diversity index through calendar years for every
site in the input compilation or dataset. By default plots Shannon `H`,
but any index supported by
[`vt_alpha_diversity()`](https://robustecologies.github.io/VerteTIME/reference/vt_alpha_diversity.md)
(or any pre-computed
[`vt_diversity()`](https://robustecologies.github.io/VerteTIME/reference/vt_diversity.md)
tibble) can be passed.

## Usage

``` r
vt_plot_alpha_timeline(
  x,
  index = "H",
  colour_by = c("taxonomic_focus", "realm")
)
```

## Arguments

- x:

  A `vt_compilation`, `vt_dataset`, or a pre-computed
  [`vt_diversity()`](https://robustecologies.github.io/VerteTIME/reference/vt_diversity.md).

- index:

  Character index name. Ignored when `x` is already a `vt_diversity`
  object.

- colour_by:

  Character; one of `"taxonomic_focus"` (default) or `"realm"` to colour
  line groups by a categorical attribute joined from the `datasets`
  table.

## Value

A `ggplot` object.

## References

Magurran, A. E., et al. (2019). *Long-term datasets reveal hidden
patterns of biodiversity change*. Nature Ecology & Evolution, 3,
1391-1394.
[doi:10.1038/s41559-019-0967-2](https://doi.org/10.1038/s41559-019-0967-2)
.

## See also

[`vt_alpha_diversity()`](https://robustecologies.github.io/VerteTIME/reference/vt_alpha_diversity.md),
[`vt_evenness()`](https://robustecologies.github.io/VerteTIME/reference/vt_evenness.md)

## Examples

``` r
if (FALSE) { # \dontrun{
co <- vt_ingest_all()
vt_plot_alpha_timeline(co, index = "H")
} # }
```
