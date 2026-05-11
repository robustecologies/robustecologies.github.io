# Construct a `vt_dataset` object

Container for a single VerteTIME dataset comprising one or more sites,
the species observed at those sites, the year-by-year abundance values,
and any environmental covariates that share the same temporal index. The
object is the unit of ingestion (one source dataset produces one
`vt_dataset`) and is the building block of
[`vt_compilation()`](https://robustecologies.github.io/VerteTIME/reference/vt_compilation.md).

## Usage

``` r
vt_dataset(
  dataset_meta,
  sites,
  species,
  observations,
  covariates = NULL,
  provenance = NULL
)

# S3 method for class 'vt_dataset'
print(x, ...)

# S3 method for class 'vt_dataset'
summary(object, ...)

# S3 method for class 'vt_dataset'
plot(
  x,
  type = c("whittaker", "rac", "alpha_timeline", "community_summary"),
  kind = NULL,
  ...
)
```

## Arguments

- dataset_meta:

  One-row data frame from the `datasets` table.

- sites:

  Data frame with one row per site for this dataset.

- species:

  Data frame with one row per species observed.

- observations:

  Long-format data frame with columns `site_id`, `species_id`, `year`,
  `value`.

- covariates:

  Optional long-format data frame with columns `site_id`, `year`,
  `covariate_name`, `value`.

- provenance:

  One-row data frame from the `data_provenance` table.

- x:

  A `vt_dataset`.

- ...:

  Forwarded to the dispatched plotter.

- object:

  A `vt_dataset`.

- type:

  One of `"whittaker"`, `"rac"`, `"alpha_timeline"`,
  `"community_summary"`. Selects the dispatched plotter.

- kind:

  Deprecated alias of `type`. Retained for back compatibility.

## Value

An object of class `vt_dataset` with the listed components stored as
`data.table`s.

For `summary`, a `tibble` with one row per species and the columns
`n_obs`, `first_year`, `last_year`, `span`, `mean_value`, `cv_value`.

A `ggplot` object (invisibly returned, drawn as side effect).

## Details

The constructor performs minimal validation; full validation is
delegated to
[`vt_validate()`](https://robustecologies.github.io/VerteTIME/reference/vt_validate.md)
for separation of concerns. The class carries `print`, `summary` and
`plot` methods. The `plot` method dispatches to one of `whittaker`,
`rac`, `alpha_timeline` or `community_summary` based on the `kind`
argument.

## References

Magurran, A. E. (2013). *Measuring Biological Diversity*. John Wiley &
Sons. ISBN 9780632056330.

## See also

[`vt_compilation()`](https://robustecologies.github.io/VerteTIME/reference/vt_compilation.md),
`print.vt_dataset()`, `summary.vt_dataset()`, `plot.vt_dataset()`,
[`vt_alpha_diversity()`](https://robustecologies.github.io/VerteTIME/reference/vt_alpha_diversity.md),
[`vt_beta_diversity()`](https://robustecologies.github.io/VerteTIME/reference/vt_beta_diversity.md),
[`vt_read()`](https://robustecologies.github.io/VerteTIME/reference/vt_read.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(vertetime)
d <- vt_read(vertetime$datasets$dataset_id[1])
print(d); summary(d); plot(d, type = "whittaker")
} # }
```
