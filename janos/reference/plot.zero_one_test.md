# Plot method for zero_one_test

Renders a scatter of the per-frequency correlations \\K(c_i)\\ against
the random frequencies \\c_i \in (\pi/5, 4\pi/5)\\ used by
[`zero_one_test`](https://robustecologies.github.io/janos/reference/zero_one_test.md).
Horizontal reference lines mark \\K = 0\\ (regular dynamics), \\K = 1\\
(chaotic dynamics) and the empirical median \\K\\ of the test.

## Usage

``` r
# S3 method for class 'zero_one_test'
plot(x, title = NULL, ...)
```

## Arguments

- x:

  A `zero_one_test` object returned by
  [`zero_one_test`](https://robustecologies.github.io/janos/reference/zero_one_test.md).

- title:

  Optional character scalar overriding the default plot title.

- ...:

  Ignored.

## Value

A `ggplot` object.

## Details

Scatter plot of the 0-1 test translation correlations

For each random frequency \\c_i\\ the 0-1 test builds the translation
variables \\(p_c, q_c)\\ and computes the correlation \\K(c_i)\\ between
the mean-square displacement and the step count. Tight clustering near
one indicates chaotic dynamics, tight clustering near zero indicates
regular dynamics, and a broad spread across the unit interval flags
either a borderline case, an undersampled trajectory or an inappropriate
stride. The scatter plot is therefore the most informative diagnostic of
the test, strictly more so than the median \\K\\ alone.

## References

Gottwald, G. A., & Melbourne, I. (2004). A new test for chaos in
deterministic systems. *Proceedings of the Royal Society A*, 460(2042),
603-611.
[doi:10.1098/rspa.2003.1183](https://doi.org/10.1098/rspa.2003.1183)

Gottwald, G. A., & Melbourne, I. (2009). On the implementation of the
0-1 test for chaos. *SIAM Journal on Applied Dynamical Systems*, 8(1),
129-145. [doi:10.1137/080718851](https://doi.org/10.1137/080718851)

## See also

[`zero_one_test`](https://robustecologies.github.io/janos/reference/zero_one_test.md),
[`print.zero_one_test`](https://robustecologies.github.io/janos/reference/print.zero_one_test.md)

## Examples

``` r
if (FALSE) { # \dontrun{
run <- dyn_sim(lorenz, t_max = 200, discard_transient = 50,
               solver = solver_rk45())
z <- zero_one_test(run)
plot(z)
} # }
```
