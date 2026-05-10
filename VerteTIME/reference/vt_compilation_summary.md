# Compilation-level summary table

Returns the `datasets` master table augmented with the count of
observation rows and the per-dataset coverage statistics.

## Usage

``` r
vt_compilation_summary(x)
```

## Arguments

- x:

  A `vt_compilation`.

## Value

A `tibble` with one row per dataset.

## References

Dornelas, M., et al. (2018). BioTIME: a database of biodiversity time
series for the Anthropocene. *Global Ecology and Biogeography*, 27(7),
760-786. [doi:10.1111/geb.12729](https://doi.org/10.1111/geb.12729) .

## See also

[`vt_compilation()`](https://robustecologies.github.io/VerteTIME/reference/vt_compilation.md),
[`vt_community_summary()`](https://robustecologies.github.io/VerteTIME/reference/vt_community_summary.md)

## Examples

``` r
if (FALSE) { # \dontrun{
co <- vt_ingest_all()
s <- vt_compilation_summary(co)
} # }
```
