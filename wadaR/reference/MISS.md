# Multiple imputation sequential sampling

This function performs multiple imputation (MI) on a numeric matrix by
sequentially sampling variables with missing values, given all other
variables in the data set.

## Usage

``` r
MISS(X, Iterations = 100, Algorithm = "GS", Fit = NULL, verbose = TRUE)
```

## Arguments

- X:

  A numeric matrix of data that contains both observed and missing
  values. Data set \\\textbf{X}\\ must not have any rows or columns that
  are completely missing. \\\textbf{X}\\ must not have any constants.
  The user must apply any data transformations appropriate for these
  models. Missing values are assumed to be Missing At Random (MAR).

- Iterations:

  The number of iterations to perform sequential sampling via MCMC
  algorithms.

- Algorithm:

  The MCMC algorithm defaults to the Gibbs Sampler (GS).

- Fit:

  An optional object of class `miss`. When supplied, `MISS` will
  continue where it left off, provided the user does not change the
  algorithm (different methods are used with different algorithms, so
  model parameters will not match). Changing algorithms requires
  starting from scratch.

- verbose:

  Logical. When `FALSE`, only the iteration prints to the console. When
  `TRUE`, which is the default, both the iteration and which variable is
  being imputed are printed to the console.

## Value

An object of class `miss` that is a list with five components:

- Algorithm:

  The algorithm that was selected.

- Imp:

  A \\M \times T\\ matrix of \\M\\ missing values and \\T\\ iterations
  that contains imputations.

- parm:

  A list of length \\J\\ for \\J\\ variables, where each component
  contains parameters associated with the prediction of missing values
  for that variable.

- PostMode:

  A vector of posterior modes. If the user intends to replace missing
  values in a data set with only one estimate per missing value (single,
  not multiple imputation), then this vector contains these values.

- Type:

  A vector of length \\J\\ for \\J\\ variables that indicates the type
  of each variable, as MISS will consider it. When `Type=1`, the
  variable is considered to be continuous. When `Type=2`, only two
  discrete values (0 and 1) were found.

## Details

Imputation is a family of statistical methods for replacing missing
values with estimates. Introduced by Rubin and Schenker (1986) and Rubin
(1987), Multiple Imputation (MI) is a family of imputation methods that
includes multiple estimates, and therefore includes variability of the
estimates.

The Multiple Imputation Sequential Sampler (MISS) function performs MI
by determining the type of variable and therefore the sampler for each
variable, and then sequentially progresses through each variable in the
data set that has missing values, updating its prediction of those
missing values given all other variables in the data set each iteration.

MI is best performed within a model, where it is called full-likelihood
imputation. Examples may be found in the "Examples" vignette. However,
sometimes it is impractical to impute within a model when there are
numerous missing values and a large number of parameters are therefore
added. As an alternative, MI may be performed on the data set before the
data is passed to the model, such as in the
[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
or
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md)
function. This is less desirable, but MISS is available for MCMC-based
MI in this case.

Missing values are initially set to column means for continuous
variables, and are set to one for discrete variables.

MISS uses the following methods and MCMC algorithms: missing values of
continuous variables are estimated with a ridge-stabilized linear
regression Gibbs sampler; missing values of binary variables that have
only 0 or 1 for values are estimated with a binary robit (t-link
logistic regression model) Gibbs sampler of Albert and Chib (1993);
missing values of discrete variables with 3 or more (ordered or
unordered) discrete values are considered continuous.

In the presence of big data, it is suggested that the user sequentially
read in batches of data that are small enough to be manageable, and then
apply the MISS function to each batch. Each batch should be
representative of the whole, of course.

It is common for multiple imputation functions to handle variable
transformations. MISS does not transform variables, but imputes what it
gets. For example, if a user has a variable that should be positive
only, then it is recommended here that the user log-transform the
variable, pass the data set to MISS, and when finished, exponentiate
both the observed and imputed values of that variable.

The `CenterScale` function should also be considered to speed up
convergence.

It is hoped that MISS is helpful, though it is not without limitation
and there are numerous alternatives outside of the `lucifer` package. If
MISS does not fulfill the needs of the user, then the following packages
are recommended: Amelia, mi, or mice. MISS emphasizes MCMC more than
these alternatives, though MISS is not as extensive. When a data set
does not have a simple structure, such as merely continuous or binary or
unordered discrete, the
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
function should be considered, where a user can easily specify
complicated structures such as multilevel, spatial or temporal
dependence, and more.

Matrix inversions are required in the Gibbs sampler. Matrix inversions
become more cumbersome as the number \\J\\ of variables increases.

If a large number of iterations is used, then the user may consider
studying the imputations for approximate convergence with the
[`BMK.Diagnostic`](https://robustecologies.github.io/lucifer/reference/BMK.Diagnostic.md)
function, by supplying the transpose of `Fit$Imp`. In the presence of
numerous missing values, say more than 100, the user may consider
iterating through the study of the imputations of 100 missing values at
a time.

## References

Albert, J.H. and Chib, S. (1993). "Bayesian Analysis of Binary and
Polychotomous Response Data". *Journal of the American Statistical
Association*, 88(422), p. 669–679.

Rubin, D.B. (1987). "Multiple Imputation for Nonresponse in Surveys".
John Wiley and Sons: New York, NY.

Rubin, D.B. and Schenker, N. (1986). "Multiple Imputation for Interval
Estimation from Simple Random Samples with Ignorable Nonresponse".
*Journal of the American Statistical Association*, 81, p. 366–374.

## See also

[`ABB`](https://robustecologies.github.io/lucifer/reference/ABB.md),
[`BMK.Diagnostic`](https://robustecologies.github.io/lucifer/reference/BMK.Diagnostic.md),
[`CenterScale`](https://robustecologies.github.io/lucifer/reference/CenterScale.md),
[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
and
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
### Create Data
N <- 20 #Number of Simulated Records
J <- 5 #Number of Simulated Variables
pM <- 0.25 #Percent Missing
Sigma <- as.positive.definite(matrix(runif(J*J),J,J))
X <- rmvn(N, rep(0,J), Sigma)
m <- sample.int(N*J, round(pM*N*J))
X[m] <- NA
head(X)

### Begin Multiple Imputation
Fit <- MISS(X, Iterations=100, Algorithm="GS", verbose=TRUE)
Fit
summary(Fit)
plot(Fit)
plot(BMK.Diagnostic(t(Fit$Imp)))

### Continue Updating if Necessary
Fit <- MISS(X, Iterations=100, Algorithm="GS", Fit, verbose=TRUE)
summary(Fit)
plot(Fit)

### Replace Missing Values in Data Set with Posterior Modes
Ximp <- X
Ximp[which(is.na(X))] <- Fit$PostMode
} # }
```
