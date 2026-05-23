# Truncated stick-breaking

The `Stick` function provides the utility of truncated stick-breaking
regarding the vector \\\theta\\. Stick-breaking is commonly referred to
as a stick-breaking process, and is used often in a Dirichlet process
(Sethuraman, 1994). It is commonly associated with infinite-dimensional
mixtures, but in practice, the "infinite" number is truncated to a
finite number, since it is impossible to estimate an infinite number of
parameters (Ishwaran and James, 2001).

## Usage

``` r
Stick(theta)
```

## Arguments

- theta:

  A vector of length \\(M-1)\\ regarding \\M\\ mixture components.

## Value

A probability vector wherein each element relates to a mixture
component.

## Details

The Dirichlet process (DP) is a stochastic process used in Bayesian
nonparametric modeling, most commonly in DP mixture models, otherwise
known as infinite mixture models. A DP is a distribution over
distributions; each draw from a DP is itself a discrete distribution. A
DP is an infinite-dimensional generalization of Dirichlet distributions.
It is called a DP because it has Dirichlet-distributed,
finite-dimensional, marginal distributions, just as the Gaussian process
has Gaussian-distributed, finite-dimensional, marginal distributions.
Distributions drawn from a DP cannot be described using a finite number
of parameters, thus the classification as a nonparametric model. The
truncated stick-breaking (TSB) process is associated with a truncated
Dirichlet process (TDP).

An example of a TSB process is cluster analysis, where the number of
clusters is unknown and treated as mixture components. In such a model,
the TSB process calculates probability vector \\\pi\\ from \\\theta\\,
given a user-specified maximum number of clusters to explore as \\C\\,
where \\C\\ is the length of \\\theta + 1\\. Vector \\\pi\\ is assigned
a TSB prior distribution (for more information, see
[`dStick`](https://robustecologies.github.io/lucifer/reference/dist.Stick.md)).

Elsewhere, each element of \\\theta\\ is constrained to the interval
(0,1), and the original TSB form is beta-distributed with the \\\alpha\\
parameter of the beta distribution constrained to 1 (Ishwaran and James,
2001). The \\\beta\\ hyperparameter in the beta distribution is usually
gamma-distributed.

A larger value for a given \\\theta_m\\ is associated with a higher
probability of the associated mixture component; however, the proportion
changes according to the position of the element in the \\\theta\\
vector.

A variety of stick-breaking processes exist. For example, rather than
each \\\theta\\ being beta-distributed, there have been other forms
introduced such as logistic and probit, among others.

## References

Ishwaran, H. and James, L. (2001). "Gibbs Sampling Methods for Stick
Breaking Priors". *Journal of the American Statistical Association*,
96(453), p. 161–173.

Sethuraman, J. (1994). "A Constructive Definition of Dirichlet Priors".
*Statistica Sinica*, 4, p. 639–650.

## See also

[`ddirichlet`](https://robustecologies.github.io/lucifer/reference/dist.Dirichlet.md),
[`dmvpolya`](https://robustecologies.github.io/lucifer/reference/dist.Multivariate.Polya.md),
and
[`dStick`](https://robustecologies.github.io/lucifer/reference/dist.Stick.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving Stick
} # }
```
