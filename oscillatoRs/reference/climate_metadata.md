# Climate index metadata

Comprehensive metadata for all 43 climate oscillation indices supported
by the package, including descriptions, physical interpretations, and
data sources.

## Usage

``` r
climate_metadata
```

## Format

A tibble with the following columns:

- index:

  Character. Index identifier.

- description:

  Character. Short description.

- long_description:

  Character. Detailed physical interpretation.

- category:

  Character. Index category.

- resolution:

  Character. Temporal resolution (monthly/daily/bimonthly).

- start_year:

  Integer. Approximate start year of data.

- source_url:

  Character. Data source URL.

## See also

[`list_indices()`](https://robustecologies.github.io/oscillatoRs/reference/list_indices.md)
for interactive index listing,
[`get_index_info()`](https://robustecologies.github.io/oscillatoRs/reference/get_index_info.md)
for single index details

## Examples

``` r
if (FALSE) { # \dontrun{
data(climate_metadata)

# View ENSO indices
climate_metadata[climate_metadata$category == "ENSO", c("index", "description")]
} # }
```
