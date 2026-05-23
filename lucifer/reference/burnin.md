# Burn-in

The `burnin` function estimates the duration of burn-in in iterations
for one or more Markov chains. "Burn-in" refers to the initial portion
of a Markov chain that is not stationary and is still affected by its
initial value.

## Usage

``` r
burnin(x, method = "BMK")
```

## Arguments

- x:

  This is a vector or matrix of posterior samples for which the number
  of burn-in iterations will be estimated.

- method:

  This argument defaults to `"BMK"`, in which case stationarity is
  estimated with the
  [`BMK.Diagnostic`](https://robustecologies.github.io/lucifer/reference/BMK.Diagnostic.md)
  function. Alternatively, the
  [`Geweke.Diagnostic`](https://robustecologies.github.io/lucifer/reference/Geweke.Diagnostic.md)
  function may be used when `method="Geweke"` or the
  [`KS.Diagnostic`](https://robustecologies.github.io/lucifer/reference/KS.Diagnostic.md)
  function may be used when `method="KS"`.

## Value

The `burnin` function returns a vector equal in length to the number of
MCMC chains in `x`, and each element indicates the maximum iteration in
burn-in.

## Details

Burn-in is a colloquial term for the initial iterations in a Markov
chain prior to its convergence to the target distribution. During
burn-in, the chain is not considered to have "forgotten" its initial
value. Burn-in is not a theoretical part of MCMC, but its use is the
norm because of the need to limit the number of posterior samples due to
computer memory. If burn-in were retained rather than discarded, then
more posterior samples would have to be retained. If a Markov chain
starts anywhere close to the center of its target distribution, then
burn-in iterations do not need to be discarded.

In the
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
function, stationarity is estimated with the
[`BMK.Diagnostic`](https://robustecologies.github.io/lucifer/reference/BMK.Diagnostic.md)
function on all thinned posterior samples of each chain, beginning at
cumulative 10% intervals relative to the total number of samples, and
the lowest number in which all chains are stationary is considered the
burn-in.

The term "burn-in" originated in electronics regarding the initial
testing of component failure at the factory to eliminate initial
failures (Geyer, 2011). Although "burn-in" has been the standard term
for decades, some are referring to these as "warm-up" iterations.

## References

Geyer, C.J. (2011). "Introduction to Markov Chain Monte Carlo". In S
Brooks, A Gelman, G Jones, and M Xiao-Li (eds.), "Handbook of Markov
Chain Monte Carlo", p. 3–48. Chapman and Hall, Boca Raton, FL.

## See also

[`BMK.Diagnostic`](https://robustecologies.github.io/lucifer/reference/BMK.Diagnostic.md),
[`deburn`](https://robustecologies.github.io/lucifer/reference/deburn.md),
[`Geweke.Diagnostic`](https://robustecologies.github.io/lucifer/reference/Geweke.Diagnostic.md),
[`KS.Diagnostic`](https://robustecologies.github.io/lucifer/reference/KS.Diagnostic.md),
and
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

## Examples

``` r
if (FALSE) { # \dontrun{
library(lucifer)
x <- rnorm(1000)
burnin(x)
} # }
```
