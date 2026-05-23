# Summary method for a zero_one_test object

Reports the 0-1 test statistics beyond the median: the empirical 95%
confidence interval of \\K\\ from the bootstrap distribution over random
test frequencies, the interquartile range, the number of frequencies in
each verdict class, the effective sample size after subsampling, the
stride selected automatically or supplied manually, and the verdict
derived from the median.

## Usage

``` r
# S3 method for class 'zero_one_test'
summary(object, ...)
```

## Arguments

- object:

  A `zero_one_test` object.

- ...:

  Unused.

## Value

A list of class `summary.zero_one_test`, invisibly.

## Details

The empirical CI is the 2.5% and 97.5% quantile of the per-frequency
correlations \\K(c_i)\\ stored in `object$K_values`; it is not a
bootstrap over the trajectory but an empirical spread across random \\c
\in (\pi/5, 4\pi/5)\\. A narrow CI centred near one supports a chaotic
verdict; a narrow CI centred near zero supports a regular verdict; a
wide spread across the unit interval indicates a borderline case or an
undersampled trajectory, as discussed in Gottwald & Melbourne (2009,
Section 2.3).

## References

Gottwald, G. A., & Melbourne, I. (2009). On the implementation of the
0-1 test for chaos. *SIAM Journal on Applied Dynamical Systems*, 8(1),
129-145. [doi:10.1137/080718851](https://doi.org/10.1137/080718851)

## See also

[`zero_one_test()`](https://robustecologies.github.io/janos/reference/zero_one_test.md)
. constructor;
[`print.zero_one_test()`](https://robustecologies.github.io/janos/reference/print.zero_one_test.md)
. compact header;
[`plot.zero_one_test()`](https://robustecologies.github.io/janos/reference/plot.zero_one_test.md)
. scatter of per-frequency correlations.

## Examples

``` r
if (FALSE) { # \dontrun{
run <- dyn_sim(lorenz, t_max = 200, discard_transient = 50)
summary(zero_one_test(run))
} # }
```
