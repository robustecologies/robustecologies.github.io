# Big data

This function enables Bayesian inference with data that is too large for
computer memory (RAM) with the simplest method: reading in batches of
data (where each batch is a section of rows), applying a function to the
batch, and combining the results.

## Usage

``` r
BigData(
  file,
  nrow,
  ncol,
  size = 1,
  Method = "add",
  CPUs = 1,
  Type = "PSOCK",
  FUN,
  ...
)
```

## Arguments

- file:

  This required argument accepts a path and filename that must refer to
  a .csv file, and that must contain only a numeric matrix without a
  header, row names, or column names.

- nrow:

  This required argument accepts a scalar integer that indicates the
  number of rows in the big data matrix.

- ncol:

  This required argument accepts a scalar integer that indicates the
  number of columns in the big data matrix.

- size:

  This argument accepts a scalar integer that specifies the number of
  rows of each batch. The last batch is not required to have the same
  number of rows as the other batches. The largest possible size, and
  therefore the fewest number of batches, should be preferred.

- Method:

  This argument accepts a scalar string, defaults to "add", and
  alternatively accepts "rbind". When `Method="rbind"`, the
  user-specified function `FUN` is applied to each batch, and results
  are combined together by rows. For example, if calculating \\\mu =
  \textbf{X}\beta\\ in, say, 10 batches, then the output column vector
  \\\mu\\ is equal to the number of rows of the big data set.

- CPUs:

  This argument accepts an integer that specifies the number of central
  processing units (CPUs) of the multicore computer or computer cluster.
  This argument defaults to `CPUs=1`, in which parallel processing does
  not occur.

- Type:

  This argument specifies the type of parallel processing to perform,
  accepting either `Type="PSOCK"` or `Type="MPI"`.

- FUN:

  This required argument accepts a user-specified function that will be
  performed on each batch. The first argument in the function must be
  the data.

- ...:

  Additional arguments are used within the user-specified function.
  Additional arguments often refer to parameters.

## Value

The `BigData` function returns output that is the result of performing a
user-specified function on batches of big data. Output is a matrix, and
may have one or more column vectors.

## Details

Big data is defined loosely here as data that is too large for computer
memory (RAM). The `BigData` function uses the split-apply-combine
strategy with a big data set. The unmanageable big data set is split
into smaller, manageable pieces (batches), a function is applied to each
batch, and results are combined.

Each iteration, the `BigData` function opens a connection to a big data
set and keeps the connection open while the `scan` function reads in
each batch of data (elsewhere, batches are often referred to chunks). A
user-specified function is applied to each batch of data, the results
are combined together, the connection is closed, and the results are
returned.

As an introductory example, suppose a statistician updates a linear
regression model, but the design matrix \\\textbf{X}\\ is too large for
computer memory. Suppose the design matrix has 100 million rows, and the
statistician specifies `size=1e6`. The statistician combines dependent
variable \\\textbf{y}\\ with design matrix \\\textbf{X}\\. Each
iteration in
[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md), or
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md),
the `BigData` function sequentially reads in one million rows of the
combined data \\\textbf{X}\\, calculates expectation vector \\\mu\\, and
finally returns the sum of the log-likelihood. The sum of the
log-likelihood is added together for all batches, and returned.

There are many limitations with this function. This function is not
fast, in the sense that the entire big data set is processed in batches,
each iteration. With iterative methods, this may perform well, albeit
slowly. There are many functions that cannot be performed on batches,
though most models in the Examples vignette may easily be updated with
big data. Large matrices of samples are unaddressed, only the data.
Although many (but not all) models may be estimated, many additional
functions in this package will not work when applied after the model has
updated. Instead, a batch or random sample of data (see the
[`read.matrix`](https://robustecologies.github.io/lucifer/reference/Matrices.md)
function for sampling from big data) should be used in the usual way, in
the `Data` argument, and the `Model` function coded in the usual way
without the `BigData` function.

Parallel processing may be performed when the user specifies `CPUs` to
be greater than one, implying that the specified number of CPUs exists
and is available. Parallelization may be performed on a multicore
computer or a computer cluster. Either a Simple Network of Workstations
(SNOW) or Message Passing Interface (MPI) is used.

There have been several alternative approaches suggested for big data.
Huang and Gelman (2005) propose that the user creates batches by
sampling from big data, updating a separate Bayesian model on each
batch, and combining the results into a consensus posterior. Scott et
al. (2013) propose a method that they call Consensus Monte Carlo.
Balakrishnan and Madigan (2006) introduced a Sequential Monte Carlo
(SMC) sampler designed for big data. Welling and Teh (2011) proposed the
stochastic gradient Langevin dynamics (SGLD) algorithm, available in the
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
function.

## References

Balakrishnan, S. and Madigan, D. (2006). "A One-Pass Sequential Monte
Carlo Method for Bayesian Analysis of Massive Datasets". *Bayesian
Analysis*, 1(2), p. 345–362.

Huang, Z. and Gelman, A. (2005) "Sampling for Bayesian Computation with
Large Datasets". *SSRN eLibrary*.

Scott, S.L., Blocker, A.W. and Bonassi, F.V. (2013). "Bayes and Big
Data: The Consensus Monte Carlo Algorithm". In *Bayes 250*.

Welling, M. and Teh, Y.W. (2011). "Bayesian Learning via Stochastic
Gradient Langevin Dynamics". *Proceedings of the 28th International
Conference on Machine Learning (ICML)*, p. 681–688.

## See also

[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`lucifer.RAM`](https://robustecologies.github.io/lucifer/reference/lucifer.RAM.md),
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md),
[`PMC.RAM`](https://robustecologies.github.io/lucifer/reference/PMC.RAM.md),
[`read.matrix`](https://robustecologies.github.io/lucifer/reference/Matrices.md),
and
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
### Below is an example of a linear regression model specification
### function in which BigData reads in a batch of 1,000 records of
### Data$N records from a data set that is too large to fully open
### in memory.

library(lucifer)
N <- 10000
J <- 10 #Number of predictors, including the intercept
X <- matrix(1,N,J)
for (j in 2:J) {X[,j] <- rnorm(N,runif(1,-3,3),runif(1,0.1,1))}
beta.orig <- runif(J,-3,3)
e <- rnorm(N,0,0.1)
y <- as.vector(tcrossprod(beta.orig, X) + e)
mon.names <- c("LP","sigma")
parm.names <- as.parm.names(list(beta=rep(0,J), log.sigma=0))
PGF <- function(Data) return(c(rnormv(Data$J,0,0.01),
     log(rhalfcauchy(1,1))))
MyData <- list(J=J, PGF=PGF, N=N, mon.names=mon.names,
     parm.names=parm.names)
filename <- tempfile("X.csv")
write.table(cbind(y,X), filename, sep=",", row.names=FALSE,
  col.names=FALSE)

Model <- function(parm, Data)
     {
     ### Parameters
     beta <- parm[1:Data$J]
     sigma <- exp(parm[Data$J+1])
     ### Log(Prior Densities)
     beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
     sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
     ### Log-Likelihood
     LL <- BigData(file=filename, nrow=Data$N, ncol=Data$J+1, size=1000,
          Method="add", CPUs=1, Type="PSOCK",
          FUN=function(x, beta, sigma) sum(dnorm(x[,1], tcrossprod(x[,-1],
               t(beta)), sigma, log=TRUE)), beta, sigma)
     ### Log-Posterior
     LP <- LL + beta.prior + sigma.prior
     Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(LP,sigma),
          yhat=0, parm=parm)
     return(Modelout)
     }

### From here, the user may update the model as usual.
} # }
```
