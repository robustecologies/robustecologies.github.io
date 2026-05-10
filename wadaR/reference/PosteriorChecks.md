# Posterior checks

Not to be confused with posterior predictive checks, this function
provides additional information about the marginal posterior
distributions of continuous parameters, such as the probability that
each posterior coefficient of the parameters (referred to generically as
\\\theta\\), is greater than zero \[\\p(\theta \> 0)\\\], the estimated
number of modes, the kurtosis and skewness of the posterior
distributions, the burn-in of each chain (for MCMC only), integrated
autocorrelation time, independent samples per minute, and acceptance
rate. A posterior correlation matrix is provided only for objects of
class `demonoid` or `pmc`.

For discrete parameters, see the
[`Hangartner.Diagnostic`](https://robustecologies.github.io/lucifer/reference/Hangartner.Diagnostic.md).

## Usage

``` r
PosteriorChecks(x, Parms = NULL)
```

## Arguments

- x:

  An object of class `demonoid`, `iterquad`, `laplace`, `pmc`, or `vb`.

- Parms:

  A vector of quoted strings to be matched for selecting parameters.
  Defaults to `NULL` and selects every parameter. Each quoted string is
  matched to one or more parameter names with the `grep` function. For
  example, if the user specifies `Parms=c("eta", "tau")`, and if the
  parameter names are beta\[1\], beta\[2\], eta\[1\], eta\[2\], and tau,
  then all parameters will be selected, because the string `eta` is
  within `beta`. Since `grep` is used, string matching uses regular
  expressions, so beware of meta-characters, though these are
  acceptable: ".", "\[", and "\]".

## Value

An object of class `posteriorchecks` that is a list with the following
components:

- Posterior.Correlation:

  A correlation matrix of the parameters selected with the `Parms`
  argument. Returned as `NA` for objects of classes `"laplace"` or
  `"vb"`.

- Posterior.Summary:

  A matrix in which each row is a parameter and there are eight columns:
  p(theta \> 0), N.Modes, Kurtosis, Skewness, Burn-In, IAT, ISM, and AR.

## Details

`PosteriorChecks` is a supplemental function that returns a list with
two components. Following is a summary of popular uses of the
`PosteriorChecks` function.

First (and only for MCMC users), the user may be considering the current
MCMC algorithm versus others. In this case, the `PosteriorChecks`
function is often used to find the two MCMC chains with the highest
[`IAT`](https://robustecologies.github.io/lucifer/reference/IAT.md), and
these chains are studied for non-randomness with a joint trace plot, via
the
[`joint.density.plot`](https://robustecologies.github.io/lucifer/reference/joint.density.plot.md)
function. The best algorithm has the chains with the highest independent
samples per minute (ISM).

Posterior correlation may be studied between model updates as well as
after a model seems to have converged. While frequentists consider
multicollinear predictor variables, Bayesians tend to consider posterior
correlation of the parameters. Models with multicollinear parameters
take more iterations to converge. Hierarchical models often have high
posterior correlations. Posterior correlation often contributes to a
lower effective sample size
([`ESS`](https://robustecologies.github.io/lucifer/reference/ESS.md)).
Common remedies include transforming the predictors, re-parameterization
to reduce posterior correlation, using WIPs (Weakly-Informative Priors),
or selecting a different numerical approximation algorithm. An example
of re-parameterization is to constrain related parameters to sum to
zero. Another approach is to specify the parameters according to a
multivariate distribution that is assisted by estimating a covariance
matrix. Some algorithms are more robust to posterior correlation than
others. Posterior correlation may be plotted with the
[`plotMatrix`](https://robustecologies.github.io/lucifer/reference/plotMatrix.md)
function, and may be useful for blocking parameters. For more
information on blockwise sampling, see the
[`Blocks`](https://robustecologies.github.io/lucifer/reference/Blocks.md)
function.

After a user is convinced of the applicability of the current MCMC
algorithm, and that the chains have converged, `PosteriorChecks` is
often used to identify multimodal marginal posterior distributions for
further study or model re-specification.

Although many marginal posterior distributions appear normally
distributed, there is no such assumption. Nonetheless, a marginal
posterior distribution tends to be distributed the same as its prior
distribution. If a parameter has a prior specified with a Laplace
distribution, then the marginal posterior distribution tends also to be
Laplace-distributed. In the common case of normality, kurtosis and
skewness may be used to identify discrepancies between the prior and
posterior, and perhaps this should be called a "prior-posterior check".

Lastly, parameter importance may be considered, in which case it is
recommended to be considered simultaneously with variable importance
from the
[`Importance`](https://robustecologies.github.io/lucifer/reference/Importance.md)
function.

## See also

[`AcceptanceRate`](https://robustecologies.github.io/lucifer/reference/AcceptanceRate.md),
[`Blocks`](https://robustecologies.github.io/lucifer/reference/Blocks.md),
[`burnin`](https://robustecologies.github.io/lucifer/reference/burnin.md),
[`ESS`](https://robustecologies.github.io/lucifer/reference/ESS.md),
[`Hangartner.Diagnostic`](https://robustecologies.github.io/lucifer/reference/Hangartner.Diagnostic.md),
[`joint.density.plot`](https://robustecologies.github.io/lucifer/reference/joint.density.plot.md),
[`IAT`](https://robustecologies.github.io/lucifer/reference/IAT.md),
[`Importance`](https://robustecologies.github.io/lucifer/reference/Importance.md),
[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`Modes`](https://robustecologies.github.io/lucifer/reference/Mode.md),
[`plotMatrix`](https://robustecologies.github.io/lucifer/reference/plotMatrix.md),
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md), and
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
### See the lucifer function for an example.
} # }
```
