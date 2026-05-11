# Compute alpha-diversity indices per site and year

Computes one or more alpha-diversity indices per
`(dataset_id, site_id, year)` combination using the long-format
`observations` table. All indices are derived in pure R with no
dependency on `vegan`, although `vegan` is available through `Suggests`
for the user's own comparisons.

## Usage

``` r
vt_alpha_diversity(x, indices = c("S", "H", "q1"))
```

## Arguments

- x:

  A `vt_compilation` or `vt_dataset`.

- indices:

  Character vector naming the indices to compute. One or more of: `"S"`
  (richness), `"H"` (Shannon entropy in nats), `"D"` (Simpson's
  inverse-concentration `1 - sum(p^2)`), `"q0"`, `"q1"`, `"q2"` (Hill
  numbers of order 0, 1, 2), `"Chao1"` (bias-corrected Chao1 estimator).
  Defaults to `c("S","H","q1")`.

## Value

A
[`vt_diversity()`](https://robustecologies.github.io/VerteTIME/reference/vt_diversity.md)
tibble keyed by `(dataset_id, site_id, year)` with columns `index` and
`value`. Hill `q1` equals `exp(H)`, Hill `q2` equals `1 / (1 - D)`, Hill
`q0` equals `S`. Chao1 returns `NA` when there are no singletons.

## Details

Hill numbers are unified diversity-of-order-q measures introduced by
Hill (1973) and recommended for cross-study comparison because they
share a common unit (effective number of species). Chao1 is the
bias-corrected estimator (Chao 1984) that recovers an estimate of the
true number of species in the assemblage from the count of singletons
and doubletons.

Site-years with `n_obs < 2` evaluate richness (`S`) but return `NA` for
entropy-based indices, since `H` is not defined on a single observation.

## References

Hill, M. O. (1973). *Diversity and evenness: a unifying notation and its
consequences*. Ecology, 54(2), 427-432.
[doi:10.2307/1934352](https://doi.org/10.2307/1934352) .

Chao, A. (1984). *Nonparametric estimation of the number of classes in a
population*. Scandinavian Journal of Statistics, 11(4), 265-270.

Magurran, A. E. (2013). *Measuring Biological Diversity*. John Wiley &
Sons. ISBN 9780632056330.

## See also

[`vt_evenness()`](https://robustecologies.github.io/VerteTIME/reference/vt_evenness.md),
[`vt_beta_diversity()`](https://robustecologies.github.io/VerteTIME/reference/vt_beta_diversity.md),
[`vt_diversity()`](https://robustecologies.github.io/VerteTIME/reference/vt_diversity.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(vertetime); co <- vertetime
a  <- vt_alpha_diversity(co, indices = c("S","H","q1","Chao1"))
summary(a)
} # }
```
