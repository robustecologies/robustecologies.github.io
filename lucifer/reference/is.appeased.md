# Appeased

Returns `TRUE` if the fitted model object passes all convergence
criteria evaluated by
[`Consort`](https://robustecologies.github.io/lucifer/reference/Consort.md),
and `FALSE` otherwise. Supports all inference families: MCMC
(`demonoid`), variational Bayes (`vb`), Laplace (`laplace`), iterative
quadrature (`iterquad`), PMC (`pmc`), and SMC (`smc`).

## Usage

``` r
is.appeased(x)
```

## Arguments

- x:

  A fitted model object of class `demonoid`, `vb`, `laplace`,
  `iterquad`, `pmc`, or `smc`.

## Value

A logical value. `TRUE` if the object passes all convergence criteria,
`FALSE` otherwise.

## Details

This is a convenience wrapper around
`Consort(x, verbose = FALSE)$appeased`. For MCMC fits, it checks the
non-adaptive requirement, acceptance rate, MCSE, ESS (bulk and tail),
Rhat, divergences, and diminishing adaptation. For non-MCMC families,
family-specific criteria are evaluated.

## See also

[`Consort`](https://robustecologies.github.io/lucifer/reference/Consort.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# MCMC fit
fit <- lucifer(Model, Data, Initial.Values,
  Iterations = 10000, Status = 1000, Thinning = 10,
  Algorithm = "CHARM", Specs = NULL)
is.appeased(fit)

# VB fit
vb_fit <- VariationalBayes(Model, Data, parm = IV,
  Iterations = 1000, Method = "Pathfinder", sir = TRUE)
is.appeased(vb_fit)
} # }
```
