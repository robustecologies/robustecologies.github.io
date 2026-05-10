# Enrich the species table with class/order/family from the taxonomy cache

Merges the GBIF backbone resolution cached at
`data-raw/_yaml/_taxonomy.tsv` (produced by `data-raw/fetch-taxonomy.R`)
into the compilation species table, populating the `class`, `order`,
`family` and `is_vertebrate` columns.

## Usage

``` r
enrich_taxonomy(spp, data_raw = here::here("data-raw"))
```

## Arguments

- spp:

  `data.table` of species with at least a `species_id` column.

- data_raw:

  Path to the data-raw root.

## Value

The species table augmented with backbone columns. Species absent from
the cache are left with `NA` in those fields.

## Details

If the taxonomy cache is missing, the function is a no-op. The
resolution is offline and uses a TSV the maintainer regenerates by
running `data-raw/fetch-taxonomy.R` after compilation changes.

## References

Chamberlain, S., et al. (2024). *rgbif: Interface to the Global
Biodiversity Information Facility API*. R package version 3.8.
<https://CRAN.R-project.org/package=rgbif>.

## See also

[`vt_ingest_all()`](https://robustecologies.github.io/VerteTIME/reference/vt_ingest_all.md),
[`vt_compilation()`](https://robustecologies.github.io/VerteTIME/reference/vt_compilation.md)

## Examples

``` r
if (FALSE) { # \dontrun{
co <- vt_ingest_all()
table(co$species$class, useNA = "ifany")
} # }
```
