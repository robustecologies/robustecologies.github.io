# Construct a `vt_compilation` object

A compilation is a collection of one or more `vt_dataset`s with a shared
schema. Constructed by binding the per-dataset relational tables
(`datasets`, `sites`, `species`, `observations`, `covariates`,
`data_provenance`) into five master tables. The compilation is the
working object for every package analytical function.

## Usage

``` r
vt_compilation(
  datasets,
  sites,
  species,
  observations,
  covariates = NULL,
  data_provenance = NULL
)

# S3 method for class 'vt_compilation'
print(x, ...)

# S3 method for class 'vt_compilation'
summary(object, ...)

# S3 method for class 'vt_compilation'
plot(x, type = c("world", "timeline", "alpha", "beta"), kind = NULL, ...)
```

## Arguments

- datasets:

  Data frame, one row per dataset.

- sites:

  Data frame, one row per site.

- species:

  Data frame, one row per species.

- observations:

  Data frame, long format, the workhorse.

- covariates:

  Optional data frame.

- data_provenance:

  Data frame, one row per dataset, the audit table.

- x:

  A `vt_compilation`.

- ...:

  Forwarded to the dispatched plotter.

- object:

  A `vt_compilation`.

- type:

  One of `"world"`, `"timeline"`, `"alpha"`, `"beta"`. Selects the
  dispatched plotter.

- kind:

  Deprecated alias of `type`. Retained for back compatibility.

## Value

An object of class `vt_compilation` with the six tables as
`data.table`s.

A `tibble` with one row per dataset and the columns of the `datasets`
table joined with `n_observation_rows`.

A `ggplot` object (invisibly returned, drawn as side effect).

## Details

All cross-table relationships use the `*_id` naming convention. The
constructor sets keys on `data.table`s for fast joins. Use
[`vt_filter()`](https://robustecologies.github.io/VerteTIME/reference/vt_filter.md)
to subset by predicate; the filter cascades referential integrity.

## References

Wilkinson, M. D., Dumontier, M., Aalbersberg, IJ. J., Appleton, G.,
Axton, M., Baak, A., Blomberg, N., Boiten, J.-W., et al. (2016). The
FAIR Guiding Principles for scientific data management and stewardship.
*Scientific Data*, 3, 160018.
[doi:10.1038/sdata.2016.18](https://doi.org/10.1038/sdata.2016.18) .

## See also

[`vt_dataset()`](https://robustecologies.github.io/VerteTIME/reference/vt_dataset.md),
[`vt_filter()`](https://robustecologies.github.io/VerteTIME/reference/vt_filter.md),
`print.vt_compilation()`, `summary.vt_compilation()`,
`plot.vt_compilation()`,
[`vt_publish()`](https://robustecologies.github.io/VerteTIME/reference/vt_publish.md),
[vertetime](https://robustecologies.github.io/VerteTIME/reference/vertetime.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(vertetime)
print(vertetime); summary(vertetime); plot(vertetime, type = "world")
} # }
```
