# Approximate Bayesian bootstrap

This function performs multiple imputation (MI) with the Approximate
Bayesian Bootstrap (ABB) of Rubin and Schenker (1986).

## Usage

``` r
ABB(X, K = 1)
```

## Arguments

- X:

  This is a vector or matrix of data that must include both observed and
  missing values. When `X` is a matrix, missing values must occur
  somewhere in the set, but are not required to occur in each variable.

- K:

  This is the number of imputations.

## Value

This function returns a list with \\K\\ components, one for each set of
imputations. Each component contains a vector of imputations equal in
length to the number of missing values in the data.

`ABB` does not currently return the mean of the imputations, or the
between-imputation variance or within-imputation variance.

## Details

The Approximate Bayesian Bootstrap (ABB) is a modified form of the
[`BayesianBootstrap`](https://robustecologies.github.io/lucifer/reference/BayesianBootstrap.md)
(Rubin, 1981) that is used for multiple imputation (MI). Imputation is a
family of statistical methods for replacing missing values with
estimates. Introduced by Rubin and Schenker (1986) and Rubin (1987), MI
is a family of imputation methods that includes multiple estimates, and
therefore includes variability of the estimates.

The data, \\\textbf{X}\\, are assumed to be independent and identically
distributed (IID), contain both observed and missing values, and its
missing values are assumed to be ignorable (meaning enough information
is available in the data that the missingness mechanism can be ignored,
if the information is used properly) and Missing Completely At Random
(MCAR). When `ABB` is used in conjunction with a propensity score
(described below), missing values may be Missing At Random (MAR).

`ABB` does not add auxiliary information, but performs imputation with
two sampling (with replacement) steps. First,
\\\textbf{X}^\star\_{obs}\\ is sampled from \\\textbf{X}\_{obs}\\. Then,
\\\textbf{X}^\star\_{mis}\\ is sampled from \\\textbf{X}^\star\_{obs}\\.
The result is a sample of the posterior predictive distribution of
\\(\textbf{X}\_{mis}\|\textbf{X}\_{obs})\\. The first sampling step is
also known as hotdeck imputation, and the second sampling step changes
the variance. Since auxiliary information is not included, `ABB` is
appropriate for missing values that are ignorable and MCAR.

Auxiliary information may be included in the process of imputation by
introducing a propensity score (Rosenbaum and Rubin, 1983; Rosenbaum and
Rubin, 1984), which is an estimate of the probability of missingness.
The propensity score is often the result of a binary logit model, where
missingness is predicted as a function of other variables. The
propensity scores are discretized into quantile-based groups, usually
quintiles. Each quintile must have both observed and missing values.
`ABB` is applied to each quintile. This is called within-class
imputation. It is assumed that the missing mechanism depends only on the
variables used to estimate the propensity score.

With \\K=1\\, `ABB` may be used in MCMC, such as in
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
more commonly along with a propensity score for missingness. MI is
performed, despite \\K=1\\, because imputation occurs at each MCMC
iteration. The practical advantage of this form of imputation is the
ease with which it may be implemented. For example, full-likelihood
imputation should perform better, but requires a chain to be updated for
each missing value.

An example of a limitation of `ABB` with propensity scores is to
consider imputing missing values of income from age in a context where
age and income have a positive relationship, and where the highest
incomes are missing systematically. `ABB` with propensity scores should
impute these highest missing incomes given the highest observed ages,
but is unable to infer beyond the observed data.

ABB has been extended (Parzen et al., 2005) to reduce bias, by
introducing a correction factor that is applied to the MI variance
estimate. This correction may be applied to output from `ABB`.

## References

Parzen, M., Lipsitz, S.R., and Fitzmaurice, G.M. (2005). "A Note on
Reducing the Bias of the Approximate Bayesian Bootstrap Imputation
Variance Estimator". *Biometrika*, 92, 4, p. 971–974.

Rosenbaum, P.R. and Rubin, D.B. (1983). "The Central Role of the
Propensity Score in Observational Studies for Causal Effects".
*Biometrika*, 70, p. 41–55.

Rosenbaum, P.R. and Rubin, D.B. (1984). "Reducing Bias in Observational
Studies Using Subclassification in the Propensity Score". *Journal of
the American Statistical Association*, 79, p. 516–524.

Rubin, D.B. (1981). "The Bayesian Bootstrap". *Annals of Statistics*, 9,
p. 130–134.

Rubin, D.B. (1987). "Multiple Imputation for Nonresponse in Surveys".
John Wiley and Sons: New York, NY.

Rubin, D.B. and Schenker, N. (1986). "Multiple Imputation for Interval
Estimation from Simple Random Samples with Ignorable Nonresponse".
*Journal of the American Statistical Association*, 81, p. 366–374.

## See also

[`BayesianBootstrap`](https://robustecologies.github.io/lucifer/reference/BayesianBootstrap.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
and
[`MISS`](https://robustecologies.github.io/lucifer/reference/MISS.md).

## Examples

``` r
if (FALSE) { # \dontrun{
library(lucifer)

### Create Data
J <- 10 #Number of variables
m <- 20 #Number of missings
N <- 50 #Number of records
mu <- runif(J, 0, 100)
sigma <- runif(J, 0, 100)
X <- matrix(0, N, J)
for (j in 1:J) X[,j] <- rnorm(N, mu[j], sigma[j])

### Create Missing Values
M1 <- rep(0, N*J)
M2 <- sample(N*J, m)
M1[M2] <- 1
M <- matrix(M1, N, J)
X <- ifelse(M == 1, NA, X)

### Approximate Bayesian Bootstrap
imp <- ABB(X, K=1)

### Replace Missing Values in X (when K=1)
X.imp <- X
X.imp[which(is.na(X.imp))] <- unlist(imp)
X.imp
} # }
```
