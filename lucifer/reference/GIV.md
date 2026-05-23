# Generate initial values

The `GIV` function generates initial values for use with the
[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md), and
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md)
functions.

## Usage

``` r
GIV(Model, Data, n = 1000, PGF = FALSE)
```

## Arguments

- Model:

  This required argument is a model specification function. For more
  information, see
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

- Data:

  This required argument is a list of data. For more information, see
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

- n:

  This is the number of attempts to generate acceptable initial values.

- PGF:

  Logical. When `TRUE`, a Parameter-Generating Function (PGF) is
  required to be in `Data`, and `GIV` will generate initial values
  according to the user-specified PGF. This argument defaults to
  `FALSE`, in which case initial values are generated randomly without
  respect to a user-specified function.

## Value

The `GIV` function returns a vector equal in length to the number of
parameters, and each element is an initial value for the associated
parameter in `Data$parm.names`. When `GIV` fails to find acceptable
initial values, each returned element is `NA`.

## Details

Initial values are required for optimization or sampling algorithms. A
user may specify initial values, or use the `GIV` function for random
generation. Initial values determined by the user may fail to produce a
finite posterior in complicated models, and the `GIV` function is here
to help.

`GIV` has several uses. First, the
[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
and
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md)
functions use `GIV` internally if unacceptable initial values are
discovered. Second, the user may use `GIV` when developing their model
specification function, `Model`, to check for potential problems. Third,
the user may prefer to randomly generate acceptable initial values.
Lastly, `GIV` is recommended when running multiple or parallel chains
with
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
(such as for later use with
[`Gelman.Diagnostic`](https://robustecologies.github.io/lucifer/reference/Gelman.Diagnostic.md))
for dispersed starting locations. For dispersed starting locations,
`GIV` should be run once for each parallel chain, and the results should
be stored per row in a matrix of initial values.

It is strongly recommended that the user specifies a
Parameter-Generating Function (PGF), and includes this function in the
list of data. Although the PGF may be specified according to the prior
distributions (possibly considered as a Prior-Generating Function), it
is often specified with a more restricted range. For example, if a user
has a model with the following prior distributions

\$\$\beta_j \sim \mathcal{N}(0, 1000), j=1,\dots,5\$\$ \$\$\sigma \sim
\mathcal{HC}(25)\$\$

then the PGF, given the prior distributions, is

`PGF <- function(Data) return(c(rnormv(5,0,1000),rhalfcauchy(1,25)))`

However, the user may not want to begin with initial values that could
be so far from zero (as determined by the variance of 1000), and may
instead prefer

`PGF <- function(Data) return(c(rnormv(5,0,10),rhalfcauchy(1,5)))`

When `PGF=FALSE`, initial values are attempted to be constrained to the
interval \\\[-100,100\]\\. This is done to prevent numeric overflows
with parameters that are exponentiated within the model specification
function. First, `GIV` passes the upper and lower bounds of this
interval to the model, and any changed parameters are noted.

At this point, it is hoped that a non-finite posterior is not found. If
found, then the remainder of the process is random and without the
previous bounds. This can be particularly problematic in the case of,
say, initial values that are the elements of a matrix that must be
positive-definite, especially with large matrices. If a random solution
is not found, then `GIV` will fail.

If the posterior is finite and `PGF=FALSE`, then initial values are
randomly generated with a normal proposal and a small variance at the
center of the returned range of each parameter. As `GIV` fails to find
acceptable initial values, the algorithm iterates toward its maximum
number of iterations, `n`. In each iteration, the variance increases for
the proposal.

Initial values are considered acceptable only when the first two
returned components of `Model` (which are `LP` and `Dev`) are finite,
and when initial values do not change through constraints, as returned
in the fifth component of the list: `parm`.

If `GIV` fails to return acceptable initial values, then it is best to
study the model specification function. When the model is complicated,
here is a suggestion. Remove the log-likelihood, `LL`, from the equation
that calculates the logarithm of the unnormalized joint posterior
density, `LP`. For example, convert `LP <- LL + beta.prior` to
`LP <- beta.prior`. Now, maximize `LP`, which is merely the set of prior
densities, with any optimization algorithm. Replace `LL`, and run the
model with initial values that are in regions of high prior density
(preferably with `PGF=TRUE`). If this fails, then the model
specification should be studied closely, because a non-finite posterior
should (especially) never be associated with regions of high prior
density.

## See also

[`as.initial.values`](https://robustecologies.github.io/lucifer/reference/as.initial.values.md),
[`Gelman.Diagnostic`](https://robustecologies.github.io/lucifer/reference/Gelman.Diagnostic.md),
[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md), and
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
data(demonsnacks)
y <- log(demonsnacks$Calories)
X <- cbind(1, as.matrix(log(demonsnacks[,c(1,4,10)]+1)))
J <- ncol(X)
for (j in 2:J) X[,j] <- CenterScale(X[,j])

mon.names <- c("LP","sigma")
parm.names <- as.parm.names(list(beta=rep(0,J), sigma=0))
pos.beta <- grep("beta", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
     beta <- rnorm(Data$J)
     sigma <- runif(1)
     return(c(beta, sigma))
     }
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

Initial.Values <- GIV(Model, MyData, PGF=TRUE)
} # }
```
