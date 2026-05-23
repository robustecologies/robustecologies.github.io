# Thin

This function reduces the number of posterior samples by retaining every
\\k\\th sample.

## Usage

``` r
Thin(x, By = 1)
```

## Arguments

- x:

  A vector or matrix of posterior samples to be thinned.

- By:

  This argument specifies that every \\k\\th posterior sample will be
  retained, and `By` defaults to 1, meaning that thinning will not
  occur, because every sample will be retained.

## Value

A thinned matrix. When `x` is a vector, the returned object is a matrix
with 1 column.

## Details

A thinned matrix of posterior samples is a matrix in which only every
\\k\\th posterior sample (or row) in the original matrix is retained.
The act of thinning posterior samples has been criticized as throwing
away information, which is correct. However, it is common practice to
thin posterior samples, usually associated with MCMC such as
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
for two reasons. First, each chain (column vector) in a matrix of
posterior samples probably has higher autocorrelation than desired,
which reduces the effective sample size (see
[`ESS`](https://robustecologies.github.io/lucifer/reference/ESS.md) for
more information). Therefore, a thinned matrix usually contains
posterior samples that are closer to independent than an un-thinned
matrix. The other reason for the popularity of thinning is that a user
may not have the random-access memory (RAM) to store large, un-thinned
matrices of posterior samples.

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
and [`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md)
automatically thin posterior samples, deviance samples, and samples of
monitored variables, according to its own user-specified argument. The
`Thin` function is made available here, should it be necessary to thin
posterior samples outside of objects of class `demonoid` or `pmc`.

## See also

[`ESS`](https://robustecologies.github.io/lucifer/reference/ESS.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
and [`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md).

## Examples

``` r
if (FALSE) { # \dontrun{
x <- matrix(runif(100), 10, 10)
Thin(x, By=2)
} # }
```
