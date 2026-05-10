# Compute partial sums of MIMFs for scale-dependent analysis

Compute cumulative reconstructions of a multivariate empirical mode
decomposition by retaining only the modes whose characteristic frequency
exceeds a sequence of cutoffs. The resulting list of high-pass
reconstructions feeds the scale-dependent dimension and persistence
computation in
[`scale_dependent_metrics`](https://robustecologies.github.io/chamaeleon/reference/scale_dependent_metrics.md)
and corresponds to equation (6) of Alberti et al. (2023).

## Usage

``` r
memd_partial_sums(memd_obj, freq_cutoffs = NULL)
```

## Arguments

- memd_obj:

  Object of class `"memd"` returned by
  [`memd`](https://robustecologies.github.io/chamaeleon/reference/memd.md).

- freq_cutoffs:

  Numeric vector. Frequency cutoffs (in Hz, on the same scale as
  `memd_obj$frequencies`) at which to evaluate the partial
  reconstruction. If `NULL` (default), the cutoffs are taken to be the
  frequencies of the individual MIMFs, yielding one reconstruction per
  mode.

## Value

A list of numeric matrices, each of dimensions `nrow(mimfs)` by
`ncol(mimfs)`, where the \\i\\-th element is the sum of all MIMFs with
characteristic frequency strictly greater than `freq_cutoffs[i]`. Two
sentinel cases override the rule: if no mode is strictly above
`freq_cutoffs[i]` the full reconstruction (sum over all modes) is
returned; if every mode is strictly above the cutoff only the
highest-frequency mode is returned. Both sentinels reproduce the
"partial sums starting from highest frequency" convention of Alberti et
al. (2023, Eq. 6) when `freq_cutoffs` is sorted in descending order, as
is the case when the default value is used.

## Details

Let \\\\c_k(t)\\\_{k=1}^{K}\\ denote the MIMFs of a multivariate signal,
ordered from highest to lowest frequency, with characteristic
frequencies \\f_1 \> f_2 \> \cdots \> f_K\\. For a cutoff \\f_c\\ that
lies strictly between two mode frequencies the partial sum is

\$\$S(f_c; t) \\=\\ \sum\_{k:\\ f_k \> f_c} c_k(t).\$\$

Each \\S(f_c)\\ is itself a multivariate signal that retains only the
high-frequency content of the decomposition. Computing EVT metrics on
\\S(f_c)\\ for a sequence of cutoffs traces a scale-resolved profile of
the dimension and persistence consumed by
[`scale_dependent_metrics`](https://robustecologies.github.io/chamaeleon/reference/scale_dependent_metrics.md),
[`chameleon_analysis`](https://robustecologies.github.io/chamaeleon/reference/chameleon_analysis.md)
and
[`chameleon_test`](https://robustecologies.github.io/chamaeleon/reference/chameleon_test.md).

**Sentinel cases.** When \\f_c \ge f_1\\ no mode is strictly above the
cutoff and the function returns the full reconstruction \\\sum_k c_k\\,
treated as the "all-modes-included" anchor of the partial-sum sequence.
When \\f_c \< f_K\\ every mode is strictly above the cutoff and the
function returns only \\c_1\\ (the highest-frequency mode), the
"minimum-scale" anchor.

## References

Alberti T et al. (2023). Chameleon attractors in turbulent flows. Chaos,
Solitons and Fractals 168:113195.
[doi:10.1016/j.chaos.2023.113195](https://doi.org/10.1016/j.chaos.2023.113195)

Rehman N, Mandic DP (2010). Multivariate empirical mode decomposition.
Proceedings of the Royal Society A 466:1291-1302.
[doi:10.1098/rspa.2009.0502](https://doi.org/10.1098/rspa.2009.0502)

Rilling G, Flandrin P, Goncalves P, Lilly JM (2007). Bivariate empirical
mode decomposition. IEEE Signal Processing Letters 14:936-939.
[doi:10.1109/LSP.2007.904710](https://doi.org/10.1109/LSP.2007.904710)

## See also

[`memd`](https://robustecologies.github.io/chamaeleon/reference/memd.md)
for the underlying decomposition,
[`scale_dependent_metrics`](https://robustecologies.github.io/chamaeleon/reference/scale_dependent_metrics.md)
for the scale-resolved EVT metrics that consume these partial
reconstructions,
[`chameleon_analysis`](https://robustecologies.github.io/chamaeleon/reference/chameleon_analysis.md)
for the full pipeline.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(42)
n <- 1000
t <- seq(0, 10, length.out = n)
x <- cbind(
  sin(2 * pi * 5 * t)   + 0.5 * sin(2 * pi * 0.5 * t),
  sin(2 * pi * 5 * t)   + 0.3 * sin(2 * pi * 0.5 * t)
) + matrix(rnorm(2 * n, sd = 0.05), ncol = 2)

decomp <- memd(x, n_directions = 16, max_imf = 6)
partial <- memd_partial_sums(decomp)
length(partial)
dim(partial[[1]])
} # }
```
