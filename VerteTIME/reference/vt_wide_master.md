# Build the wide-format master series table for the compilation

Produces an LPD-style wide table where each row is one observation
series (one `(dataset_id, site_id, species_id)` triple), the leading
columns are all the metadata anchored to that triple (dataset metadata,
site metadata, species metadata, primary-reference fields), and the
trailing columns are the calendar years from `year_min` to `year_max`
populated with the abundance values (`NA` where unobserved). The result
is the canonical analytic-ready denormalisation of the compilation and a
drop-in for any tool that expects the LPD wide layout.

## Usage

``` r
vt_wide_master(x, year_min = NULL, year_max = NULL)
```

## Arguments

- x:

  A `vt_compilation`.

- year_min:

  Optional integer; defaults to the compilation-wide minimum.

- year_max:

  Optional integer; defaults to the compilation-wide maximum.

## Value

A `data.table` with one row per series and `length(year_min:year_max)`
trailing year columns.

## Details

The leading metadata block carries `dataset_id`, `site_id`,
`species_id`, `genus`, `species_epithet`, `class`, `order`, `family`,
`is_vertebrate`, the dataset-level columns (`taxonomic_focus`, `realm`,
`biome`, `ecoregion`, `country_iso3`, `unit_class`,
`is_community_metric_eligible`, `n_species`, `n_years`), the site-level
columns (`latitude_dd`, `longitude_dd`, `coord_precision_km`, `habitat`,
`elevation_m`), and the provenance columns
(`primary_reference_citation`, `primary_reference_doi`,
`primary_reference_kind`).

The compilation is materialised in the canonical `VT_NNN` namespace; the
wide master inherits the same ids without further transformation.

## References

McGill, B. J., et al. (2007). *Species abundance distributions: moving
beyond single prediction theories to integration within an ecological
framework*. Ecology Letters, 10(10), 995-1015.
[doi:10.1111/j.1461-0248.2007.01094.x](https://doi.org/10.1111/j.1461-0248.2007.01094.x)
.

## See also

[`vt_publish()`](https://robustecologies.github.io/VerteTIME/reference/vt_publish.md),
[`vt_export()`](https://robustecologies.github.io/VerteTIME/reference/vt_export.md),
[`vt_compilation()`](https://robustecologies.github.io/VerteTIME/reference/vt_compilation.md)

## Examples

``` r
if (FALSE) { # \dontrun{
co <- vt_ingest_all()
wm <- vt_wide_master(co)
dim(wm)                 # ~556 rows; ~150 columns
names(wm)[1:25]
} # }
```
