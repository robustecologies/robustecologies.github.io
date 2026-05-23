# Model specification time

The `Model.Spec.Time` function returns the time in minutes to evaluate a
model specification a number of times, as well as the evaluations per
minute, and componentwise iterations per minute.

## Usage

``` r
Model.Spec.Time(Model, Initial.Values, Data, n = 1000)
```

## Arguments

- Model:

  A model specification function. For more information, see
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

- Initial.Values:

  A vector of initial values for the parameters.

- Data:

  A list of data. For more information, see
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

- n:

  The number of evaluations of the model specification, and accuracy
  increases with `n`.

## Value

A list with three components:

- Time:

  The time in minutes to evaluate the model specification `n` times.

- Evals.per.Minute:

  The number of evaluations of the model specification per minute: `n` /
  `Time`. This is also the number of iterations per minute in an
  algorithm that is not componentwise, where one evaluation occurs per
  iteration.

- Componentwise.Iters.per.Minute:

  The number of iterations per minute in a componentwise MCMC algorithm,
  such as AMWG or MWG. It is the evaluations per minute divided by the
  number of parameters, since an evaluation must occur for each
  parameter, for each iteration.

## Details

The largest single factor to affect the run-time of an algorithm,
whether it is with
[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md), or
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md),
is the time that it takes to evaluate the model specification. This has
also been observed in Rosenthal (2007). It is highly recommended that a
user of the `lucifer` package should attempt to reduce the run-time of
the model specification, usually by testing alternate forms of code for
speed. This is especially true with big data, such as with the
[`BigData`](https://robustecologies.github.io/lucifer/reference/BigData.md)
function.

Every function in the lucifer package is byte-compiled, which is a
recent option in R. This reduces run-time, thanks to Tierney's compiler
package in base R. The model specification, however, is specified by the
user, and should be byte-compiled. The reduction in run-time may range
from mild to dramatic, depending on the model. It is highly recommended
that users concerned with run-time activate the compiler package and use
the `cmpfun` function, as per the example below.

A model specification function that is optimized for speed and involves
many records may result in a model update run-time that is considerably
less than other popular MCMC-based software algorithms that loop through
records, even when those algorithms are coded in `C` or other fast
languages. For a comparison, see the "Lucifer Tutorial" vignette.

However, if a model specification function in the lucifer package is not
fully vectorized (contains `for` loops and `apply` functions), then
run-time will typically be slower than other popular MCMC-based software
algorithms.

The speed of calculating the model specification function is affected by
parameter constraints, such as with the
[`interval`](https://robustecologies.github.io/lucifer/reference/interval.md)
function. Parameter constraints may add considerable variability to the
speed of this calculation, and usually more variation occurs with
initial values that are far from the target distributions.

Distributions in the `lucifer` package usually have logical checks to
ensure correctness. These checks may slow the calculation of the
density, for example. If the user is confident these checks are
unnecessary for their model, then the user may copy the code to a new
function name and comment-out the checks to improve speed.

When speed is of paramount importance, a user is encouraged to
experiment with writing the model specification function in another
language such as in `C++` with the `Rcpp` package, and calling drop-in
replacement functions from within the `Model` function, or re-writing
the model function entirely in `C++`.

Another use for `Model.Spec.Time` is to allow the user to make an
informed decision about which MCMC algorithm to select, given the speed
of their model specification. For example, the Adaptive
Metropolis-within-Gibbs (AMWG) of Roberts and Rosenthal (2009) is
currently the adaptive MCMC algorithm of choice in general in many
cases, but this choice is conditional on run-time. While other MCMC
algorithms in `lucifer` evaluate the model specification function once
per iteration, componentwise algorithms such as in the MWG family
evaluate it once per parameter per iteration, significantly increasing
run-time per iteration in large models. The `Model.Spec.Time` function
may forewarn the user of the associated run-time, and if it should be
decided to go with an alternate MCMC algorithm, such as AMM, in which
each element of its covariance matrix must stabilize for the chains to
become stationary. AMM, for example, will require many more iterations
of burn-in (for more information, see the
[`burnin`](https://robustecologies.github.io/lucifer/reference/burnin.md)
function), but with numerous iterations, allows more thinning. A general
recommendation may be to use AMWG when `Componentwise.Iters.per.Minute`
\>= 1000, but this is subjective and best determined by each user for
each model. A better decision may be made by comparing inference methods
with the
[`Arena`](https://robustecologies.github.io/lucifer/reference/Arena.md)
function for a particular model.

## References

Roberts, G.O. and Rosenthal, J.S. (2009). "Examples of Adaptive MCMC".
*Computational Statistics and Data Analysis*, 18, p. 349–367.

## See also

[`BigData`](https://robustecologies.github.io/lucifer/reference/BigData.md),
[`interval`](https://robustecologies.github.io/lucifer/reference/interval.md),
[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`Arena`](https://robustecologies.github.io/lucifer/reference/Arena.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md), and
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
library(lucifer)
data(demonsnacks)
y <- log(demonsnacks$Calories)
X <- cbind(1, as.matrix(log(demonsnacks[,c(1,4,10)]+1)))
J <- ncol(X)
for (j in 2:J) {X[,j] <- CenterScale(X[,j])}

mon.names <- "LP"
parm.names <- as.parm.names(list(beta=rep(0,J), sigma=0))
pos.beta <- grep("beta", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) return(c(rnormv(Data$J,0,10), rhalfcauchy(1,5)))
MyData <- list(J=J, PGF=PGF, X=X, mon.names=mon.names,
     parm.names=parm.names, pos.beta=pos.beta, pos.sigma=pos.sigma, y=y)

Model <- function(parm, Data)
     {
     beta <- parm[Data$pos.beta]
     sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
     parm[Data$pos.sigma] <- sigma
     beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
     sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
     mu <- tcrossprod(Data$X, t(beta))
     LL <- sum(dnorm(Data$y, mu, sigma, log=TRUE))
     LP <- LL + beta.prior + sigma.prior
     Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
          yhat=rnorm(length(mu), mu, sigma), parm=parm)
     return(Modelout)
     }

set.seed(666)
Initial.Values <- GIV(Model, MyData, PGF=TRUE)
Model.Spec.Time(Model, Initial.Values, MyData)
} # }
```
