# Pareto smoothed importance sampling

Applies Pareto smoothed importance sampling (PSIS) to a matrix of
log-importance-ratios, stabilizing the upper tail via GPD fitting. This
standalone function is useful outside the LOO-CV context, for example
when combining multiple importance sampling distributions.

## Usage

``` r
PSIS(log_ratios)
```

## Arguments

- log_ratios:

  An \\N \times S\\ matrix of log-importance-ratios.

## Value

An object of class `psis`, which is a list with components:

- log_weights:

  An \\N \times S\\ matrix of normalized smoothed log-weights (rows sum
  to 0 on the log scale, i.e. weights sum to 1).

- pareto_k:

  Numeric vector of Pareto shape parameters.

## Details

Implementation of `PSIS`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

## References

Vehtari, A., Simpson, D., Gelman, A., Yao, Y., and Gabry, J. (2024).
"Pareto smoothed importance sampling." *Journal of Machine Learning
Research*, 25(72), 1–58.

## See also

[`LOO`](https://robustecologies.github.io/lucifer/reference/LOO.md)

## Examples

``` r
if (FALSE) { # \dontrun{
log_ratios <- matrix(rnorm(100 * 500), nrow = 100, ncol = 500)
result <- PSIS(log_ratios)
str(result)
} # }
```
