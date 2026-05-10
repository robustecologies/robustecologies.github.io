# Rejection sampling

The `RejectionSampling` function implements rejection sampling of a
target density given a proposal density.

## Usage

``` r
RejectionSampling(
  Model,
  Data,
  mu,
  S,
  df = Inf,
  logc,
  n = 1000,
  CPUs = 1,
  Type = "PSOCK"
)
```

## Arguments

- Model:

  A model specification function. For more information, see
  [`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md).

- Data:

  A list of data. For more information, see
  [`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md).

- mu:

  A mean vector \\\mu\\ for the multivariate normal or multivariate t
  proposal density.

- S:

  A covariance matrix \\\Sigma\\ for the multivariate normal or
  multivariate t proposal density.

- df:

  A scalar degrees of freedom parameter \\\nu\\. It defaults to
  infinity, in which case the multivariate normal density is used.

- logc:

  The logarithm of the rejection sampling constant.

- n:

  The number of independent draws to be simulated from the proposal
  density.

- CPUs:

  An integer that specifies the number of central processing units
  (CPUs) of the multicore computer or computer cluster. This argument
  defaults to `CPUs=1`, in which parallel processing does not occur.

- Type:

  The type of parallel processing to perform, accepting either
  `Type="PSOCK"` or `Type="MPI"`.

## Value

The `RejectionSampling` function returns an object of class `rejection`,
which is a matrix of accepted, independent, simulated draws from the
target distribution.

## Details

Rejection sampling (von Neumann, 1951) is a Monte Carlo method for
drawing independent samples from a distribution that is proportional to
the target distribution, \\f(x)\\, given a sampling distribution,
\\g(x)\\, from which samples can readily be drawn, and for which there
is a finite constant \\c\\.

Here, the target distribution, \\f(x)\\, is the result of the `Model`
function. The sampling distribution, \\g(x)\\, is either a multivariate
normal or multivariate t-distribution. The parameters of \\g(x)\\ (`mu`,
`S`, and `df`) are used to create random draws, \\\theta\\, of the
sampling distribution, \\g(x)\\. These draws, \\\theta\\, are used to
evaluate the target distribution, \\f(x)\\, via the `Model`
specification function. The evaluations of the target distribution,
sampling distribution, and the constant are used to create a probability
of acceptance for each draw, by comparing to a vector of \\n\\ uniform
draws, \\u\\. Each draw, \\\theta\\ is accepted if \$\$u \le
\frac{f(\theta\|\textbf{y})}{cg(\theta)}\$\$

Before beginning rejection sampling, a goal of the user is to find the
bounding constant, \\c\\, such that \\f(\theta\|\textbf{y}) \le
cg(\theta)\\ for all \\\theta\\. These are all expressed in logarithms,
so the goal is to find \\\log f(\theta\|\textbf{y}) - \log g(\theta) \le
\log(c)\\ for all \\\theta\\. This is done by maximizing \\\log
f(\theta\|\textbf{y}) - \log g(\theta)\\ over all \\\theta\\. By using,
say,
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md)
to find the modes of the parameters of interest, and using the resultant
`LP`, the mode of the logarithm of the joint posterior distribution, as
\\\log(c)\\.

The `RejectionSampling` function performs one iteration of rejection
sampling. Rejection sampling is often iterated, then called the
rejection sampling algorithm, until a sufficient number or proportion of
\\\theta\\ is accepted. An efficient rejection sampling algorithm has a
high acceptance rate. However, rejection sampling becomes less efficient
as the model dimension (the number of parameters) increases.

Extensions of rejection sampling include Adaptive Rejection Sampling
(ARS) (either derivative-based or derivative-free) and Adaptive
Rejection Metropolis Sampling (ARMS), as in Gilks et al. (1995). The
random-walk Metropolis algorithm (Metropolis et al., 1953) combined the
rejection sampling (a method of Monte Carlo simulation) of von Neumann
(1951) with Markov chains.

Parallel processing may be performed when the user specifies `CPUs` to
be greater than one, implying that the specified number of CPUs exists
and is available. Parallelization may be performed on a multicore
computer or a computer cluster. Either a Simple Network of Workstations
(SNOW) or Message Passing Interface (MPI) is used. With small data sets
and few samples, parallel processing may be slower, due to computer
network communication. With larger data sets and more samples, the user
should experience a faster run-time.

This function is similar to the `rejectsampling` function in the
`LearnBayes` package.

## References

Gilks, W.R., Best, N.G., Tan, K.K.C. (1995). "Adaptive Rejection
Metropolis Sampling within Gibbs Sampling". Journal of the Royal
Statistical Society. Series C (Applied Statistics), Vol. 44, No. 4, p.
455–472.

Metropolis, N., Rosenbluth, A.W., Rosenbluth, M.N., and Teller, E.
(1953). "Equation of State Calculations by Fast Computing Machines".
Journal of Chemical Physics, 21, p. 1087-1092.

von Neumann, J. (1951). "Various Techniques Used in Connection with
Random Digits. Monte Carlo Methods". National Bureau Standards, 12, p.
36–38.

## See also

[`dmvn`](https://robustecologies.github.io/lucifer/reference/dist.Multivariate.Normal.md),
[`dmvt`](https://robustecologies.github.io/lucifer/reference/dist.Multivariate.t.md),
[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
and
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## Suppose an output object of class laplace is called Fit:
rs <- RejectionSampling(Model, MyData, mu=Fit$Summary1[,1],
     S=Fit$Covar, df=Inf, logc=Fit$LP.Final, n=1000)
rs
} # }
```
