# Compute evenness indices per site and year

Returns Pielou's `J = H / log(S)` and (optionally) Simpson evenness
`E_D = (1/sum(p^2)) / S`.

## Usage

``` r
vt_evenness(x, kind = c("pielou", "simpson"))
```

## Arguments

- x:

  A `vt_compilation` or `vt_dataset`.

- kind:

  Character, one of `"pielou"` (default) or `"simpson"`.

## Value

A
[`vt_diversity()`](https://robustecologies.github.io/VerteTIME/reference/vt_diversity.md)
tibble keyed by `(dataset_id, site_id, year)` with columns `index` and
`value`.

## Details

Evenness is undefined when richness is zero or one; the function returns
`NA` in that case rather than `0` or `Inf`.

## References

Pielou, E. C. (1966). *The measurement of diversity in different types
of biological collections*. Journal of Theoretical Biology, 13, 131-144.
[doi:10.1016/0022-5193(66)90013-0](https://doi.org/10.1016/0022-5193%2866%2990013-0)
.

## See also

[`vt_alpha_diversity()`](https://robustecologies.github.io/VerteTIME/reference/vt_alpha_diversity.md),
[`vt_diversity()`](https://robustecologies.github.io/VerteTIME/reference/vt_diversity.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(vertetime); co <- vertetime
e  <- vt_evenness(co, kind = "pielou")
summary(e)
} # }
```
