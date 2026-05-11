# VerteTIME single-dataset demo compilation

A minimal one-dataset slice of
[vertetime](https://robustecologies.github.io/VerteTIME/reference/vertetime.md)
(the first canonical identifier, `VT_001`). Used by vignette and example
chunks where loading the full compilation would be wasteful.

## Usage

``` r
vt_demo
```

## Format

A `vt_compilation` carrying the same six slots as
[vertetime](https://robustecologies.github.io/VerteTIME/reference/vertetime.md),
restricted to a single dataset.

## Source

Slice of
[vertetime](https://robustecologies.github.io/VerteTIME/reference/vertetime.md)
produced at build time. Provenance flows through
`vt_demo$data_provenance` and reaches the same primary-literature
citations as the parent compilation.

## References

Wilkinson, M. D., Dumontier, M., Aalbersberg, IJ. J., Appleton, G.,
Axton, M., Baak, A., Blomberg, N., Boiten, J.-W., et al. (2016). *The
FAIR Guiding Principles for scientific data management and stewardship*.
Scientific Data, 3, 160018.
[doi:10.1038/sdata.2016.18](https://doi.org/10.1038/sdata.2016.18) .

## See also

[vertetime](https://robustecologies.github.io/VerteTIME/reference/vertetime.md),
[`vt_compilation()`](https://robustecologies.github.io/VerteTIME/reference/vt_compilation.md),
[`vt_read()`](https://robustecologies.github.io/VerteTIME/reference/vt_read.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(vt_demo)
print(vt_demo); summary(vt_demo)
plot(vt_demo, type = "world")
} # }
```
