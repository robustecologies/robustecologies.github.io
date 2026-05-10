# Kullback-Leibler divergence

This function calculates the Kullback-Leibler divergence (KLD) between
two probability distributions, with applications in lowest posterior
loss probability intervals, posterior predictive checks, prior
elicitation, reference priors, and Variational Bayes.

## Usage

``` r
KLD(px, py, base = exp(1))
```

## Arguments

- px:

  a numeric vector of probability densities (or log-densities)
  representing \\p(\textbf{x})\\.

- py:

  a numeric vector of probability densities (or log-densities)
  representing \\p(\textbf{y})\\. If log-densities are used, both `px`
  and `py` must be log-densities.

- base:

  the logarithmic base. Defaults to `exp(1)` (natural units, nats); use
  `base = 2` for binary units (bits).

## Value

A list with the following components:

- KLD.px.py:

  elementwise \\KLD_i(p(x_i) \|\| p(y_i))\\

- KLD.py.px:

  elementwise \\KLD_i(p(y_i) \|\| p(x_i))\\

- mean.KLD:

  the elementwise mean of the two above

- sum.KLD.px.py:

  directed divergence \\KLD(p(x) \|\| p(y))\\

- sum.KLD.py.px:

  directed divergence \\KLD(p(y) \|\| p(x))\\

- mean.sum.KLD:

  mean of the two directed divergences

- intrinsic.discrepancy:

  minimum of the two directed divergences

## Details

KLD is an asymmetric measure of the divergence between two probability
distributions \\p(\textbf{y})\\ and \\p(\textbf{x})\\ (Kullback and
Leibler, 1951). Here, \\p(\textbf{y})\\ represents the true or reference
distribution, and \\p(\textbf{x})\\ represents a model or approximation.
In Bayesian inference, KLD measures information gain in moving from
prior to posterior, and is the basis of reference priors and lowest
posterior loss intervals
([`LPL.interval`](https://robustecologies.github.io/lucifer/reference/LPL.interval.md)).

## References

Kullback, S. and Leibler, R.A. (1951). "On Information and Sufficiency".
*The Annals of Mathematical Statistics*, 22(1), p. 79–86.

Bernardo, J.M. (2005). "Intrinsic Credible Regions: An Objective
Bayesian Approach to Interval Estimation". *Sociedad de Estadistica e
Investigacion Operativa*, 14(2), p. 317–384.

## See also

[`LPL.interval`](https://robustecologies.github.io/lucifer/reference/LPL.interval.md)
and
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(lucifer)
px <- dnorm(runif(100), 0, 1)
py <- dnorm(runif(100), 0.1, 0.9)
KLD(px, py)
} # }
```
