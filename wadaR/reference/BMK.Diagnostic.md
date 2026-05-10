# BMK convergence diagnostic

Given a matrix of posterior samples from MCMC, the `BMK.Diagnostic`
function calculates Hellinger distances between consecutive batches for
each chain. This is useful for monitoring convergence of MCMC chains.

## Usage

``` r
BMK.Diagnostic(X, batches = 10)
```

## Arguments

- X:

  a matrix of posterior samples or an object of class `demonoid`, in
  which case it uses the posterior samples in `X$Posterior1`.

- batches:

  the number of batches on which the convergence diagnostic will be
  calculated. Defaults to 10.

## Value

An object of class `bmk` that is a \\J \times B\\ matrix of Hellinger
distances between consecutive batches for \\J\\ parameters of posterior
samples. The number of columns, \\B\\, is equal to the number of batches
minus one. The `BMK.Diagnostic` function is similar to the `bmkconverge`
function in package BMK.

## Details

Hellinger distance is used to quantify dissimilarity between two
probability distributions. It is based on the Hellinger integral,
introduced by Hellinger (1909). Traditionally, Hellinger distance is
bound to the interval \[0,1\], though another popular form occurs in the
interval \[0,\\\sqrt{2}\\\]. A higher value of Hellinger distance is
associated with more dissimilarity between the distributions.

Convergence is assumed when Hellinger distances are below a threshold,
indicating that posterior samples are similar between consecutive
batches. If all Hellinger distances beyond a given batch of samples is
below the threshold, then `burnin` is suggested to occur immediately
before the first batch of satisfactory Hellinger distances.

As an aid to interpretation, consider a matrix of 1,000 posterior
samples from three chains: `beta[1]`, `beta[2]`, and `beta[3]`. With 10
batches, the column names are: 100, 200, ..., 900. A Hellinger distance
for the chain `beta[1]` at 100 is the Hellinger distance between two
batches: samples 1-100, and samples 101:200.

A benefit to using `BMK.Diagnostic` is that the resulting Hellinger
distances may easily be plotted with the `plotMatrix` function, allowing
the user to see quickly which consecutive batches of which chains were
dissimilar. This makes it easier to find problematic chains.

The `BMK.Diagnostic` is calculated automatically in the
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
function, and is one of the criteria in the
[`Consort`](https://robustecologies.github.io/lucifer/reference/Consort.md)
function regarding the recommendation of when to stop updating the
Markov chain Monte Carlo (MCMC) sampler in
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

For more information on the related topics of burn-in and stationarity,
see the `burnin` and
[`is.stationary`](https://robustecologies.github.io/lucifer/reference/is.stationary.md)
functions, and the accompanying vignettes.

## References

Boone, E.L., Merrick, J.R. and Krachey, M.J. (2013). "A Hellinger
Distance Approach to MCMC Diagnostics". *Journal of Statistical
Computation and Simulation*, in press.

Hellinger, E. (1909). "Neue Begrundung der Theorie quadratischer Formen
von unendlichvielen Veranderlichen" (in German). *Journal fur die reine
und angewandte Mathematik*, 136, p. 210–271.

## See also

[`burnin`](https://robustecologies.github.io/lucifer/reference/burnin.md),
[`Consort`](https://robustecologies.github.io/lucifer/reference/Consort.md),
[`is.stationary`](https://robustecologies.github.io/lucifer/reference/is.stationary.md),
and
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

## Examples

``` r
if (FALSE) { # \dontrun{
N <- 1000
J <- 10
Theta <- matrix(runif(N*J), N, J)
colnames(Theta) <- paste("beta[", 1:J, "]", sep="")
for (i in 2:N) {Theta[i,1] <- Theta[i-1,1] + rnorm(1)}
HD <- BMK.Diagnostic(Theta, batches=10)
plot(HD, title="Hellinger distance between batches")
} # }
```
