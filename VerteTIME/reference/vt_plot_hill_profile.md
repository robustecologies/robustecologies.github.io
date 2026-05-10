# Hill-number profile per dataset across orders 0 to 3

For every dataset, computes the compilation-wide Hill number for q in
the sequence `seq(0, 3, by = 0.25)` averaged across years and sites, and
plots one line per dataset coloured by `taxonomic_focus`. Hill profiles
characterise the diversity-of-order continuum: q=0 is richness, q=1 is
exp(Shannon), q=2 is inverse Simpson, large q saturates to dominance
reciprocal.

## Usage

``` r
vt_plot_hill_profile(x)
```

## Arguments

- x:

  A `vt_compilation`.

## Value

A `ggplot`.

## References

Hill (1973); Chao et al. (2014).

## See also

[`vt_alpha_diversity()`](https://robustecologies.github.io/VerteTIME/reference/vt_alpha_diversity.md),
[`vt_plot_alpha_hist()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_alpha_hist.md)

## Examples

``` r
if (FALSE)  vt_plot_hill_profile(vt_ingest_all())  # \dontrun{}
```
