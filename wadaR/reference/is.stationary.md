# Logical check of stationarity

This function returns `TRUE` if the object is stationary according to
the
[`Geweke.Diagnostic`](https://robustecologies.github.io/lucifer/reference/Geweke.Diagnostic.md)
function, and `FALSE` otherwise.

## Usage

``` r
is.stationary(x)
```

## Arguments

- x:

  A vector, matrix, or object of class `demonoid`.

## Value

A logical value indicating whether or not the supplied object is
stationary according to the
[`Geweke.Diagnostic`](https://robustecologies.github.io/lucifer/reference/Geweke.Diagnostic.md)
function.

## Details

Stationarity, here, refers to the limiting distribution in a Markov
chain. A series of samples from a Markov chain, in which each sample is
the result of an iteration of a Markov chain Monte Carlo (MCMC)
algorithm, is analyzed for stationarity, meaning whether or not the
samples trend or its moments change across iterations. A stationary
posterior distribution is an equilibrium distribution, and assessing
stationarity is an important diagnostic toward inferring Markov chain
convergence.

In the cases of a matrix or an object of class `demonoid`, all Markov
chains (as column vectors) must be stationary for `is.stationary` to
return `TRUE`.

Alternative ways to assess stationarity of chains are to use the
[`BMK.Diagnostic`](https://robustecologies.github.io/lucifer/reference/BMK.Diagnostic.md)
or
[`Heidelberger.Diagnostic`](https://robustecologies.github.io/lucifer/reference/Heidelberger.Diagnostic.md)
functions.

## See also

[`BMK.Diagnostic`](https://robustecologies.github.io/lucifer/reference/BMK.Diagnostic.md),
[`Geweke.Diagnostic`](https://robustecologies.github.io/lucifer/reference/Geweke.Diagnostic.md),
[`Heidelberger.Diagnostic`](https://robustecologies.github.io/lucifer/reference/Heidelberger.Diagnostic.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)

## Examples

``` r
if (FALSE) { # \dontrun{
is.stationary(rnorm(100))
is.stationary(matrix(rnorm(100), 10, 10))
} # }
```
