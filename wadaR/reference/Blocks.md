# Blocks

The `Blocks` function returns a list of \\N\\ blocks of parameters, for
use with some MCMC algorithms in the
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
function. Blocks may be created either sequentially, or from a
hierarchical clustering of the posterior correlation matrix.

## Usage

``` r
Blocks(Initial.Values, N, PostCor = NULL)
```

## Arguments

- Initial.Values:

  This required argument is a vector of initial values.

- N:

  This optional argument indicates the desired number of blocks. If
  omitted, then the truncated square root of the number of initial
  values is used. If a posterior correlation matrix is supplied to
  `PostCor`, then `N` may be a scalar, or have length two. If `N` has
  length two, then the first element indicates the minimum number of
  blocks, and the second element indicates the maximum number of blocks,
  and the number of blocks is the maximum of the mean silhouette width
  for each hierarchical cluster solution.

- PostCor:

  This optional argument defaults to `NULL`, in which case sequential
  blocking is performed. If a posterior correlation matrix is supplied,
  then blocks are created based on hierarchical clustering.

## Value

The `Blocks` function returns an object of class `blocks`, which is a
list. Each component of the list is a block of parameters, and
parameters are indicated by their position in the initial values vector.

## Details

Usually, there is more than one target distribution in MCMC, in which
case it must be determined whether it is best to sample from target
distributions individually, in groups, or all at once. Blockwise
sampling (also called block updating) refers to splitting a multivariate
vector into groups called blocks, and each block is sampled separately.
A block may contain one or more parameters.

Parameters are usually grouped into blocks such that parameters within a
block are as correlated as possible, and parameters between blocks are
as independent as possible. This strategy retains as much of the
parameter correlation as possible for blockwise sampling, as opposed to
componentwise sampling where parameter correlation is ignored. The
[`PosteriorChecks`](https://robustecologies.github.io/lucifer/reference/PosteriorChecks.md)
function can be used on the output of previous runs to find highly
correlated parameters.

Advantages of blockwise sampling are that a different MCMC algorithm may
be used for each block (or parameter, for that matter), creating a more
specialized approach (though different algorithms by block are not
supported here), the acceptance of a newly proposed state is likely to
be higher than sampling from all target distributions at once in high
dimensions, and large proposal covariance matrices can be reduced in
size, which is most helpful again in high dimensions.

Disadvantages of blockwise sampling are that correlations probably exist
between parameters between blocks, and each block is updated while
holding the other blocks constant, ignoring these correlations of
parameters between blocks. Without simultaneously taking everything into
account, the algorithm may converge slowly or never arrive at the proper
solution.

Large-dimensional models often require blockwise sampling. For example,
with thousands of parameters, a componentwise algorithm must evaluate
the model specification function once per parameter per iteration,
resulting in an algorithm that may take longer than is acceptable to
produce samples. The most practical solution is to group parameters into
\\N\\ blocks, and each iteration the algorithm evaluates the model
specification function \\N\\ times, each with a reduced set of
parameters. Blockwise sampling is offered in the following algorithms:
Adaptive-Mixture Metropolis (AMM), Adaptive Metropolis-within-Gibbs
(AMWG), Automated Factor Slice Sampler (AFSS), Elliptical Slice Sampler
(ESS), Hit-And-Run Metropolis (HARM), Metropolis-within-Gibbs (MWG),
Random-Walk Metropolis (RWM), Robust Adaptive Metropolis (RAM), Slice
Sampler (Slice), and the Univariate Eigenvector Slice Sampler (UESS).

The `Blocks` function performs either a sequential assignment of
parameters to blocks when posterior correlation is not supplied, or uses
hierarchical clustering to create blocks based on posterior correlation.
If posterior correlation is supplied, then the user may specify a range
of the number of blocks to consider, and the optimal number of blocks is
considered to be the maximum of the mean silhouette width of each
hierarchical clustering. Silhouette width is calculated as per the
`cluster` package. Hierarchical clustering is performed on the distance
matrix calculated from the dissimilarity matrix (1 - abs(PostCor)) of
the posterior correlation matrix.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
and
[`PosteriorChecks`](https://robustecologies.github.io/lucifer/reference/PosteriorChecks.md).

## Examples

``` r
if (FALSE) { # \dontrun{
library(lucifer)

### Create the default number of sequentially assigned blocks:
Initial.Values <- rep(0,1000)
MyBlocks <- Blocks(Initial.Values)
MyBlocks

### Or, a pre-specified number of sequentially assigned blocks:
#Initial.Values <- rep(0,1000)
#MyBlocks <- Blocks(Initial.Values, N=20)

### If scaled diagonal covariance matrices are desired:
#VarCov <- list()
#for (i in 1:length(MyBlocks))
#  VarCov[[i]] <- diag(length(MyBlocks[[i]]))*2.38^2/length(MyBlocks[[i]])

### Or, determine the number of blocks in the range of 2 to 50 from
### hierarchical clustering on the posterior correlation matrix of an
### object, say called Fit, output from lucifer:
#MyBlocks <- Blocks(Initial.Values, N=c(2,50),
#  PostCor=cor(Fit$Posterior1))
#lapply(MyBlocks, length) #See the number of parameters per block
} # }
```
