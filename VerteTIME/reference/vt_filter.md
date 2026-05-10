# Filter a compilation by a tidy-eval predicate against the `datasets` table

Applies a tidy-evaluation predicate (the same syntax as
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html))
to the `datasets` master table and cascades the resulting `dataset_id`
set through `sites`, `observations`, `covariates` and `data_provenance`
so that the returned object remains referentially consistent.

## Usage

``` r
vt_filter(x, ...)
```

## Arguments

- x:

  A `vt_compilation`.

- ...:

  Tidy-evaluation expressions; reference any column of the `datasets`
  table (for example `taxonomic_focus == "birds"`, `n_species >= 5`,
  `year_max >= 2010`). Multiple expressions are combined with logical
  AND.

## Value

A `vt_compilation` with the filtered tables; FK integrity is preserved.

## Details

Where the user wants to filter on attributes that live in the `sites`
table (for example `latitude_dd > 0`), the recommended pattern is to
first join `sites` into a derived view, identify the relevant
`dataset_id`s, and pass those as a `dataset_id %in% ...` predicate. This
keeps the cascade simple and predictable.

## References

Wickham, H., et al. (2024). *dplyr: A Grammar of Data Manipulation*. R
package. <https://dplyr.tidyverse.org>.

## See also

[`vt_compilation()`](https://robustecologies.github.io/VerteTIME/reference/vt_compilation.md),
[`vt_export()`](https://robustecologies.github.io/VerteTIME/reference/vt_export.md),
[`vt_publish()`](https://robustecologies.github.io/VerteTIME/reference/vt_publish.md)

## Examples

``` r
if (FALSE) { # \dontrun{
co <- vt_ingest_all()
birds <- vt_filter(co, taxonomic_focus == "birds")
big   <- vt_filter(co, n_species >= 5, n_years >= 10)
} # }
```
