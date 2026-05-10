# Parameter names

This function creates a vector of parameter names from a list of
parameters, and the list may contain any combination of scalars,
vectors, matrices, upper-triangular matrices, and arrays.

## Usage

``` r
as.parm.names(x, uppertri = NULL)
```

## Arguments

- x:

  This required argument is a list of named parameters. The list may
  contain scalars, vectors, matrices, and arrays. The value of the named
  parameters does not matter here, though they are usually set to zero.
  However, if a missing value occurs, then the associated element is
  omitted in the output.

- uppertri:

  This optional argument must be a vector with a length equal to the
  number of named parameters. Each element in `uppertri` must be either
  a 0 or 1, where a 1 indicates that an upper triangular matrix will be
  used for the associated element in the vector of named parameters.
  Each element of `uppertri` is associated with a named parameter. The
  `uppertri` argument does not function with arrays.

## Value

This function returns a vector of parameter names.

## Details

Each `model` function for
[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md), or
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md)
requires a vector of parameters (specified at first as `Initial.Values`)
and a list of data. One component in the list of data must be named
`parm.names`. Each element of `parm.names` is a name associated with the
corresponding parameter in `Initial.Values`.

The `parm.names` vector is easy to program explicitly for a simple
model, but can require considerably more programming effort for more
complicated models. The `as.parm.names` function is a utility function
designed to minimize programming by the user.

For example, a simple model may only require
`parm.names <- c("alpha", "beta[1]", "beta[2]", "sigma")`. A more
complicated model may contain hundreds of parameters that are a
combination of scalars, vectors, matrices, upper-triangular matrices,
and arrays, and is the reason for the `as.parm.names` function. The code
for the above is `as.parm.names(list(alpha=0, beta=rep(0,2), sigma=0))`.

In the case of an upper-triangular matrix, simply pass the full matrix
to `as.parm.names` and indicate that only the upper-triangular will be
used via the `uppertri` argument. For example,
`as.parm.names(list(beta=rep(0,J),U=diag(K)), uppertri=c(0,1))` creates
parameter names for a vector of \\\beta\\ parameters of length \\J\\ and
an upper-triangular matrix \\\textbf{U}\\ of dimension \\K\\.

Numerous examples may be found in the accompanying "Examples" vignette.

## See also

[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md), and
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
library(lucifer)
N <- 100
J <- 5
y <- rnorm(N,0,1)
X <- matrix(runif(N*J,-2,2),N,J)
S <- diag(J)
T <- diag(2)
mon.names <- c("LP","sigma")
parm.names <- as.parm.names(list(log.sigma=0, beta=rep(0,J), S=diag(J),
     T=diag(2)), uppertri=c(0,0,0,1))
MyData <- list(J=J, N=N, S=S, T=T, X=X, mon.names=mon.names,
     parm.names=parm.names, y=y)
MyData
} # }
```
