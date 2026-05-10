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

## See also

[vertetime](https://robustecologies.github.io/VerteTIME/reference/vertetime.md),
[`vt_compilation()`](https://robustecologies.github.io/VerteTIME/reference/vt_compilation.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(vt_demo)
summary(vt_demo)
} # }
```
