# Construct a `vt_provenance` object

Wrapper of the `data_provenance` audit table with a print method that
summarises by `primary_reference_kind` and lists the per-dataset DOI
when the user asks for a single dataset.

## Usage

``` r
vt_provenance(x)

# S3 method for class 'vt_provenance'
print(x, ...)
```

## Arguments

- x:

  A `vt_provenance`.

- ...:

  Ignored.

## Value

Object of class `vt_provenance`.

## References

Wilkinson, M. D., et al. (2016). *The FAIR Guiding Principles for
scientific data management and stewardship*. Scientific Data, 3, 160018.
[doi:10.1038/sdata.2016.18](https://doi.org/10.1038/sdata.2016.18) .

## See also

[`vt_compilation()`](https://robustecologies.github.io/VerteTIME/reference/vt_compilation.md),
[`vt_publish()`](https://robustecologies.github.io/VerteTIME/reference/vt_publish.md),
[vertetime](https://robustecologies.github.io/VerteTIME/reference/vertetime.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(vertetime)
p <- vt_provenance(vertetime$data_provenance)
print(p)
} # }
```
