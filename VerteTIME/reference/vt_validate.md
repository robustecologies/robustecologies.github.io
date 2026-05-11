# Validate a `vt_dataset` or `vt_compilation`

Runs structural and semantic checks against the relational schema and
reports every failure as one row of a `tibble`. Empty output means a
clean bill of health.

## Usage

``` r
vt_validate(x, level = c("dataset", "compilation"))
```

## Arguments

- x:

  A `vt_dataset` or `vt_compilation`.

- level:

  Either `"dataset"` (default) for per-dataset checks, or
  `"compilation"` for cross-dataset uniqueness checks (DOI uniqueness,
  no `dataset_id` collisions, etc.).

## Value

A `tibble` with one row per failed check and the columns `level`,
`dataset_id`, `code`, `message`. An empty `tibble` indicates no failure.

## Details

Checks performed at the `dataset` level:

- column 1 is `year`, integer, in `[vt.year_min, vt.year_max]`;

- `species_id` matches the species regex;

- `value` is non-negative or `NA`;

- no duplicated `(site_id, species_id, year)` triples;

- `latitude_dd` in `[-90, 90]` and `longitude_dd` in `[-180, 180]`, with
  `(0, 0)` flagged as suspicious unless explicitly allowed.

Checks at the `compilation` level:

- unique `dataset_id`;

- unique `primary_reference_doi` (warning only, since some datasets
  legitimately share a primary reference when one paper reports several
  sites);

- cross-dataset `species_id` consistency (one row per `species_id`).

## References

Borer, E. T., et al. (2009). *Some simple guidelines for effective data
management*. Bulletin of the Ecological Society of America, 90(2),
205-214.
[doi:10.1890/0012-9623-90.2.205](https://doi.org/10.1890/0012-9623-90.2.205)
.

## See also

[`vt_ingest_all()`](https://robustecologies.github.io/VerteTIME/reference/vt_ingest_all.md),
[`vt_register_dataset()`](https://robustecologies.github.io/VerteTIME/reference/vt_register_dataset.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(vertetime); co <- vertetime
v  <- vt_validate(co, level = "compilation")
if (nrow(v)) print(v) else message("Compilation validation passed.")
} # }
```
