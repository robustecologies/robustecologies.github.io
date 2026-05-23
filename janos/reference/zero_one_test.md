# 0-1 test for chaos

Implements the Gottwald-Melbourne (2004, 2009) test that returns a
median correlation coefficient \\K \in \[0, 1\]\\: values close to zero
indicate regular (periodic or quasi-periodic) dynamics, values close to
one indicate chaos. The test operates on a scalar observable of a
trajectory and is independent of the system dimension and of the
embedding details.

## Usage

``` r
zero_one_test(
  x,
  var = NULL,
  n_c = 100L,
  n_cut = NULL,
  stride = NULL,
  seed = 42
)
```

## Arguments

- x:

  A `dyn_sim` object.

- var:

  Name of the scalar observable (defaults to the first state).

- n_c:

  Number of random frequencies \\c\\ to sample (default 100).

- n_cut:

  Cut-off index for the mean-square displacement curve, must be much
  smaller than `length(signal)`. Default: one tenth of the signal
  length.

- stride:

  Optional positive integer to downsample the observable before
  computing the test. Oversampled continuous flows must be subsampled so
  that consecutive points are only weakly correlated (Gottwald and
  Melbourne, 2009). When `NULL` the stride is selected automatically as
  the first lag at which \\\|\rho\_\phi(\tau)\| \< 0.5\\, capped so that
  at least 500 points survive subsampling.

- seed:

  Random seed.

## Value

An S3 object of class `zero_one_test` with components `K` (scalar, the
median), `K_values` (numeric vector), `c_values` (numeric vector),
`n_cut`, `var`.

## Details

Gottwald-Melbourne 0-1 test for chaos

For a scalar series \\\phi_n\\, the test constructs translation
variables \\p_c(n) = \sum\_{j=1}^{n} \phi_j \cos(j c)\\ and \\q_c(n) =
\sum\_{j=1}^{n} \phi_j \sin(j c)\\ at a random frequency \\c \in (0,
\pi)\\. The mean-square displacement \\M_c(n) =
\mathrm{mean}\_j\[(p_c(j+n) - p_c(j))^2 + (q_c(j+n) - q_c(j))^2\]\\
scales linearly in \\n\\ for chaotic dynamics and bounded in \\n\\ for
regular dynamics. The correlation coefficient \\K_c\\ between \\M_c(n)\\
and \\n\\ is reported; the statistic is the median \\K\\ over many
random \\c\\ values.

## References

Gottwald, G. A., & Melbourne, I. (2004). A new test for chaos in
deterministic systems. *Proceedings of the Royal Society A*, 460(2042),
603-611.
[doi:10.1098/rspa.2003.1183](https://doi.org/10.1098/rspa.2003.1183)

Gottwald, G. A., & Melbourne, I. (2009). On the implementation of the
0-1 test for chaos. *SIAM Journal on Applied Dynamical Systems*, 8(1),
129-145. [doi:10.1137/080718851](https://doi.org/10.1137/080718851)

## See also

[`print.zero_one_test`](https://robustecologies.github.io/janos/reference/print.zero_one_test.md),
[`plot.zero_one_test`](https://robustecologies.github.io/janos/reference/plot.zero_one_test.md),
[`lyapunov_spectrum`](https://robustecologies.github.io/janos/reference/lyapunov_spectrum.md),
[`correlation_dimension`](https://robustecologies.github.io/janos/reference/correlation_dimension.md),
[`poincare_section`](https://robustecologies.github.io/janos/reference/poincare_section.md),
[`bifurcation_diagram`](https://robustecologies.github.io/janos/reference/bifurcation_diagram.md)

## Examples

``` r
if (FALSE) { # \dontrun{
run <- dyn_sim(lorenz, t_max = 200, discard_transient = 50,
               solver = solver_rk45())
z <- zero_one_test(run)
print(z); plot(z)
} # }
```
