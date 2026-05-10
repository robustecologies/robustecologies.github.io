# Positional timeline of biodiversity time-series compilations

Gantt-style figure that places VerteTIME alongside the major
biodiversity time-series compilations cited in the field on the same
calendar axis. The compilations shown have heterogeneous taxonomic
scope, ranging from vertebrate-only (LPD, RivFishTIME, BBS, PECBMS, CBC,
eBird Trends, TEAM/Wildlife Insights) to multi-taxon (GPDD, BioTIME,
PREDICTS), and the figure does not claim a uniform vertebrate scope
across the rows. The `taxonomic_scope` column of
[`vt_compilations()`](https://robustecologies.github.io/VerteTIME/reference/vt_compilations.md)
disambiguates each entry. The figure is the centrepiece of section 3 of
the manuscript and the visual hook of the package homepage.

## Usage

``` r
vt_plot_database_timeline(x = NULL, extra = NULL)
```

## Arguments

- x:

  Optional `vt_compilation`; ignored except to drive the position of the
  `VerteTIME` row to the compilation's actual `(year_min, year_max)`.

- extra:

  Optional `tibble` with the same shape as
  [`vt_compilations()`](https://robustecologies.github.io/VerteTIME/reference/vt_compilations.md);
  appended to the built-in list before plotting.

## Value

A `ggplot` object.

## Details

Each row is a horizontal segment from `first_release_year` to
`latest_release_year`, coloured by `community_data` (yes / no).
VerteTIME is highlighted as the bottom row with a contrasting colour.
The caption is explicit that this figure is *positional only*: the
compilations shown are not the source of any VerteTIME series; partial
overlaps with individual VerteTIME series, where they exist, are noted
in the `data_provenance` table and are decoupled from this comparator.

## References

Halpern, B. S., et al. (2020). *Recent pace of change in human impact on
the world's ocean*. Scientific Reports, 9, 11609.
[doi:10.1038/s41598-019-47201-9](https://doi.org/10.1038/s41598-019-47201-9)
.

## See also

[`vt_compilations()`](https://robustecologies.github.io/VerteTIME/reference/vt_compilations.md),
[`vt_plot_world_map()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_world_map.md)

## Examples

``` r
if (FALSE) { # \dontrun{
p <- vt_plot_database_timeline()
print(p)
} # }
```
