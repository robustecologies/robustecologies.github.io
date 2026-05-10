# Logical check of propriety

This function provides a logical check of the propriety of a univariate
prior probability distribution or the joint posterior distribution.

## Usage

``` r
is.proper(f, a, b, tol = 1e-05)
```

## Arguments

- f:

  Either a probability density function or an object of class
  `demonoid`, `laplace`, `pmc`, or `vb`.

- a:

  The lower limit of integration, and may be negative infinity.

- b:

  The upper limit of integration, and may be positive infinity.

- tol:

  The tolerance, indicating the allowable difference from one.

## Value

A logical value indicating whether or not the univariate prior or joint
posterior probability distribution integrates to one within its
specified limits. `TRUE` is returned for a proper univariate probability
distribution.

## Details

A proper probability distribution is a probability distribution that
integrates to one, and an improper probability distribution does not
integrate to one. If a probability distribution integrates to any
positive and finite value other than one, then it is an improper
distribution, but is merely unnormalized. An unnormalized distribution
may be multiplied by a constant so that it integrates to one.

In Bayesian inference, the posterior probability distribution should be
proper. An improper prior distribution can cause an improper posterior
distribution. When the posterior distribution is improper, inferences
are invalid, it is non-integrable, and Bayes factors cannot be used
(though there are exceptions).

To avoid these problems, it is suggested that the prior probability
distribution should be proper, though it is possible to use an improper
prior distribution and have it result in a proper posterior
distribution.

To check the propriety of a univariate prior probability distribution,
create a function `f`. For example, to check the propriety of a vague
normal distribution, such as

\$\$\theta \sim \mathcal{N}(0,1000)\$\$

the function is `function(x){dnormv(x,0,1000)}`. Next, set the lower and
upper limits of integration, `a` and `b`. Internally, this function
calls `integrate` from base R, which uses adaptive quadrature. By using
\\f(x)\\ as shorthand for the specified function, `is.proper` will check
to see if the area of the following integral is one:

\$\$\int^b_a f(x)dx\$\$

Multivariate prior probability distributions currently cannot be checked
for approximate propriety. This is currently unavailable in this
package.

To check the propriety of the joint posterior distribution, the only
argument to be supplied is an object of class `demonoid`, `iterquad`,
`laplace`, `pmc`, or `vb`. The `is.proper` function checks the logarithm
of the marginal likelihood (see
[`LML`](https://robustecologies.github.io/lucifer/reference/LML.md)) for
a finite value, and returns `TRUE` when the LML is finite. This
indicates that the marginal likelihood is finite for all observed
\\\textbf{y}\\ in the model data set. This implies:

\$\$\int p(\theta\|\textbf{y})p(\theta)d\theta \< \infty\$\$

If the object is of class `demonoid` and the algorithm was adaptive, or
if the object is of class `iterquad`, `laplace`, or `vb` and the
algorithm did not converge, then `is.proper` will return `FALSE` because
LML was not estimated. In this case, it is possible for the joint
posterior to be proper, but `is.proper` will be unable to determine
propriety without the estimate of LML. If desired, the
[`LML`](https://robustecologies.github.io/lucifer/reference/LML.md) may
be estimated by the user, and if it is finite, then the joint posterior
distribution is proper.

## See also

[`dnormv`](https://robustecologies.github.io/lucifer/reference/dist.Normal.Variance.md),
[`integrate`](https://rdrr.io/r/stats/integrate.html),
[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LML`](https://robustecologies.github.io/lucifer/reference/LML.md),
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Prior probability distribution
is.proper(function(x) {dnormv(x, 0, 1000)}, -Inf, Inf) # x ~ N(0,1000)
is.proper(function(x) {dhalfcauchy(x, 25)}, 0, Inf) # x ~ HC(25)
is.proper(function(x) {dunif(x, 0, 1)}, 0, 1) # x ~ U(0,1)
is.proper(function(x) {dunif(x, -Inf, Inf)}, -Inf, Inf) # x ~ U(-Inf,Inf)

# Joint posterior distribution
# Assumes Fit is an object of class demonoid, iterquad, laplace, or pmc
# is.proper(Fit)
} # }
```
