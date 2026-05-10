# Print method for zero_one_test

Prints a compact summary of a `zero_one_test` object: the observable
used, the effective sample size and stride, the number of random test
frequencies, the median correlation \\K\\, and a verdict string derived
from \\K\\.

## Usage

``` r
# S3 method for class 'zero_one_test'
print(x, ...)
```

## Arguments

- x:

  A `zero_one_test` object returned by
  [`zero_one_test`](https://robustecologies.github.io/janos/reference/zero_one_test.md).

- ...:

  Ignored.

## Value

Invisibly returns `x`.

## Details

Print a compact summary of a zero_one_test object

The verdict classifies the test output as `"regular"` when \\K \< 0.3\\,
`"chaotic"` when \\K \> 0.7\\, and `"borderline"` otherwise. The two
cutoffs follow the practical recommendation of Gottwald and Melbourne
(2009) for subsampled flows; borderline verdicts typically signal either
insufficient data, an inappropriate stride, or quasiperiodic dynamics
that the scalar \\K\\ cannot cleanly separate.

## References

Gottwald, G. A., & Melbourne, I. (2009). On the implementation of the
0-1 test for chaos. *SIAM Journal on Applied Dynamical Systems*, 8(1),
129-145. [doi:10.1137/080718851](https://doi.org/10.1137/080718851)

## See also

[`zero_one_test()`](https://robustecologies.github.io/janos/reference/zero_one_test.md)
— constructor;
[`summary.zero_one_test()`](https://robustecologies.github.io/janos/reference/summary.zero_one_test.md)
— detailed statistics and CI;
[`plot.zero_one_test()`](https://robustecologies.github.io/janos/reference/plot.zero_one_test.md)
— scatter of per-frequency correlations;
[`lyapunov_spectrum()`](https://robustecologies.github.io/janos/reference/lyapunov_spectrum.md)
— complementary chaos diagnostic;
[`correlation_dimension()`](https://robustecologies.github.io/janos/reference/correlation_dimension.md)
— fractal-dimension diagnostic.

## Examples

``` r
if (FALSE) { # \dontrun{
run <- dyn_sim(lorenz, t_max = 200, discard_transient = 50,
               solver = solver_rk45())
z <- zero_one_test(run)
print(z)
} # }
```
