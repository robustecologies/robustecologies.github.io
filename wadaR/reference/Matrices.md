# Matrix utility functions

These are utility functions for working with matrices.

The `as.indicator.matrix` function creates an indicator matrix from a
vector. This function is useful for converting a discrete vector into a
matrix in which each column represents one of the discrete values, and
each occurence of that value in the related column is indicated by a
one, and is otherwise filled with zeroes. This function is similar to
the `class.ind` function in the nnet package.

The `as.inverse` function returns the matrix inverse of `x`. The `solve`
function in base R also returns the matrix inverse, but `solve` can
return a matrix that is not symmetric, and can fail due to
singularities. The `as.inverse` function tries to use the `solve`
function to return a matrix inverse, and when it fails due to a
singularity, `as.inverse` uses eigenvalue decomposition (in which
eigenvalues below a tolerance are replaced with the tolerance), and
coerces the result to a symmetric matrix. This is similar to the
`solvcov` function in the fpc package.

The `as.parm.matrix` function prepares a correlation, covariance, or
precision matrix in two important ways. First, `as.parm.matrix` obtains
the parameters for the matrix specified in the `x` argument by matching
the name of the matrix in the `x` argument with any parameters in
`parm`, given the parameter names in the `Data` listed in `parm.names`.
These obtained parameters are organized into a matrix as the elements of
the upper-triangular, including the diagonal. A copy is made, without
the diagonal, and the lower-triangular is filled in, completing the
matrix. Second, `as.parm.matrix` checks for positive-definiteness. If
matrix `x` is positive-definite, then the matrix is stored as a variable
called `luciferMatrix` in a new environment called `LDEnv`. If matrix
`x` is not positive-definite, then `luciferMatrix` in `LDEnv` is sought
as a replacement. If this variable exists, then it is used to replace
the matrix. If not, then the matrix is replaced with an identity matrix.
Back in the model specification, after using `as.parm.matrix`, it is
recommended that the user also pass the resulting matrix back into the
`parm` vector, so the sampler or algorithm knows that the elements of
the matrix have changed.

The `as.positive.definite` function returns the nearest
positive-definite matrix for a matrix that is square and symmetric
(Higham, 2002). This version is intended only for covariance and
precision matrices, and has been optimized for speed. A more extensible
function is `nearPD` in the matrixcalc package, which is also able to
work with correlation matrices, and matrices that are asymmetric.

The `as.positive.semidefinite` function iteratively seeks to return a
square, symmetric matrix that is at least positive-semidefinite, by
replacing each negative eigenvalue and calculating its projection. This
is intended only for covariance and precision matrices. A similar
function is `makePsd` in the RTAQ package, though it is not iterative,
and returns matrices that fail a logical check with
`is.positive.semidefinite`.

The `as.symmetric.matrix` function accepts either a vector or matrix,
and returns a symmetric matrix. In the case of a vector, it can be
either all elements of the matrix, or the lower triangular. In the case
of a `x` being entered as a matrix, this function tolerates non-finite
values in one triangle (say, the lower), as long as the corresponding
element is finite in the other (say, the upper) triangle.

The `Cov2Cor` function converts a covariance matrix into a correlation
matrix, and accepts the covariance matrix either in matrix or vector
form. This function may be useful inside a model specification and also
with converting posterior draws of the elements of a covariance matrix
to a correlation matrix. `Cov2Cor` is an expanded form of the `cov2cor`
function in the `stats` package, where `Cov2Cor` is also able to accept
and return a vectorized matrix.

The `CovEstim` function estimates a covariance matrix with one of
several methods. This is mainly used by
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
where the `parm` argument receives the posterior modes. See the `CovEst`
argument for more details.

The `GaussHermiteCubeRule` function returns a matrix of nodes and a
vector of weights for a `dims`-dimensional integral given \\N\\
univariate nodes. The number of multivariate nodes will differ from the
number of univariate nodes. This function is for use with multivariate
quadrature, often called cubature. This has been adapted from the
`multiquad` function in the NominalLogisticBiplot package. The
[`GaussHermiteQuadRule`](https://robustecologies.github.io/lucifer/reference/Math.md)
function is a univariate version. A customized univariate `rule` may be
supplied when constraints necessitate that one or more nodes and weights
had to be altered.

The `Hessian` returns a symmetric, Hessian matrix, which is a matrix of
second partial derivatives. The estimation of the Hessian matrix is
approximated numerically using Richardson extrapolation by default. This
is a slow function. This function is not intended to be called by the
user, but is made available here. This is essentially the `hessian`
function from the numDeriv package, adapted to Lucifer.

The `is.positive.definite` function is a logical test of whether or not
a matrix is positive-definite. A \\k \times k\\ symmetric matrix
\\\textbf{X}\\ is positive-definite if all of its eigenvalues are
positive (\\\lambda_i \> 0, i \in k\\). All main-diagonal elements must
be positive. The determinant of a positive-definite matrix is always
positive, so a positive-definite matrix is always nonsingular.
Non-symmetric, positive-definite matrices exist, but are not considered
here.

The `is.positive.semidefinite` function is a logical test of whether or
not a matrix is positive-semidefinite. A \\k x k\\ symmetric matrix
\\\textbf{X}\\ is positive-semidefinite if all of its eigenvalues are
non-negative (\\\lambda_i \ge 0, i \in k\\).

The `is.square.matrix` function is a logical test of whether or not a
matrix is square. A square matrix is a matrix with the same number of
rows and columns, and is usually represented as a \\k \times k\\ matrix
\\\textbf{X}\\.

The `is.symmetric.matrix` function is a logical test of whether or not a
matrix is symmetric. A symmetric matrix is a square matrix that is equal
to its transpose, \\\textbf{X} = \textbf{X}^T\\. For example, where
\\i\\ indexes rows and \\j\\ indexes columns, \\\textbf{X}\_{i,j} =
\textbf{X}\_{j,i}\\. This differs from the `isSymmetric` function in
base R that is inexact, using `all.equal`.

The `Jacobian` function estimates the Jacobian matrix, which is a matrix
of all first-order partial derivatives of the `Model`. The Jacobian
matrix is estimated by default with forward finite-differencing, or
optionally with Richardson extrapolation. This function is not intended
to be called by the user, but is made available here. This is
essentially the `jacobian` function from the numDeriv package, adapted
to lucifer.

The `logdet` function returns the logarithm of the determinant of a
positive-definite matrix via the Cholesky decomposition. The determinant
is a value associated with a square matrix, and was used historically to
*determine* if a system of linear equations has a unique solution. The
term *determinant* was introduced by Gauss, where Laplace referred to it
as the resultant. When the determinant is zero, the matrix is singular
and non-invertible; there are either no solutions or many solutions. A
unique solution exists when the determinant is non-zero. The `det`
function in base R works well for small matrices, but can return
erroneously return zero in larger matrices. It is better to work with
the log-determinant.

The `lower.triangle` function returns a vector of the lower triangular
elements of a matrix, and the diagonal is included when `diag=TRUE`.

The `read.matrix` function is provided here as one of many convenient
ways to read a numeric matrix into R. The most common method of storing
data in R is the data frame, because it is versatile. For example, a
data frame may contain character, factor, and numeric variables
together. For iterative estimation, common in Bayesian inference, the
data frame is much slower than the numeric matrix. For this reason, the
lucifer package does not use data frames, and has not traditionally
accepted character or factor data. The `read.matrix` function returns
either an entire numeric matrix, or row-wise samples from a numeric
matrix. Samples may be taken from a matrix that is too large for
available computer memory (RAM), such as with big data.

The `SparseGrid` function returns a sparse grid for a \\J\\-dimensional
integral with accuracy \\K\\, given Gauss-Hermite quadrature rules. A
grid of order \\K\\ provides an exact result for a polynomial of total
order of \\2K - 1\\ or less. `SparseGrid` returns a matrix of nodes and
a vector of weights. A sparse grid is more efficient than the full grid
in the `GaussHermiteCubeRule` function. This has been adapted from the
SparseGrid package.

The `TransitionMatrix` function has several uses. A user may supply a
vector of marginal posterior samples of a discrete Markov chain as
`theta.y`, and an observed posterior transition matrix is returned.
Otherwise, a user may supply data (`y.theta`) and/or a prior
(`p.theta`), in which case a posterior transition matrix is returned. A
common row-wise prior is the dirichlet distribution. Transition
probabilities are from row element to column element.

The `tr` function returns the trace of a matrix. The trace of a matrix
is the sum of the elements in the main diagonal of a square matrix. For
example, the trace of a \\k \times k\\ matrix \\\textbf{X}\\, is
\\\sum\_{k=1} \textbf{X}\_{k,k}\\.

The `upper.triangle` function returns a vector of the upper triangular
elements of a matrix, and the diagonal is included when `diag=TRUE`.

S3 method: apply [`as()`](https://rdrr.io/r/methods/as.html) to objects
of class `inverse`.

S3 method: apply [`as()`](https://rdrr.io/r/methods/as.html) to objects
of class `parm.matrix`.

S3 method: apply [`as()`](https://rdrr.io/r/methods/as.html) to objects
of class `positive.definite`.

S3 method: apply [`as()`](https://rdrr.io/r/methods/as.html) to objects
of class `positive.semidefinite`.

S3 method: apply [`as()`](https://rdrr.io/r/methods/as.html) to objects
of class `symmetric.matrix`.

See Details.

See Details.

See Details.

See Details.

S3 method: apply [`is()`](https://rdrr.io/r/methods/is.html) to objects
of class `positive.definite`.

S3 method: apply [`is()`](https://rdrr.io/r/methods/is.html) to objects
of class `positive.semidefinite`.

S3 method: apply [`is()`](https://rdrr.io/r/methods/is.html) to objects
of class `square.matrix`.

S3 method: apply [`is()`](https://rdrr.io/r/methods/is.html) to objects
of class `symmetric.matrix`.

See Details.

See Details.

S3 method: apply `lower()` to objects of class `triangle`.

S3 method: apply `read()` to objects of class `matrix`.

See Details.

See Details.

See Details.

S3 method: apply `upper()` to objects of class `triangle`.

## Usage

``` r
as.indicator.matrix(x)

as.inverse(x)

as.parm.matrix(
  x,
  k,
  parm,
  Data,
  a = -Inf,
  b = Inf,
  restrict = FALSE,
  chol = FALSE
)

as.positive.definite(x)

as.positive.semidefinite(x)

as.symmetric.matrix(x, k = NULL)

Cov2Cor(Sigma)

CovEstim(Model, parm, Data, Method = "Hessian")

GaussHermiteCubeRule(N, dims, rule)

Hessian(Model, parm, Data, Interval = 1e-04, Method = "Richardson")

is.positive.definite(x)

is.positive.semidefinite(x)

is.square.matrix(x)

is.symmetric.matrix(x)

Jacobian(Model, parm, Data, Interval = 1e-06, Method = "simple")

logdet(x)

lower.triangle(x, diag = FALSE)

read.matrix(
  file,
  header = FALSE,
  sep = ",",
  nrow = 0,
  samples = 0,
  size = 0,
  na.rm = FALSE
)

SparseGrid(J, K)

TransitionMatrix(theta.y = NULL, y.theta = NULL, p.theta = NULL)

tr(x)

upper.triangle(x, diag = FALSE)
```

## Arguments

- x:

  This is a matrix (though `as.symmetric.matrix` also accepts vectors).

- k:

  For `as.parm.matrix`, this is a required argument, indicating the
  dimension of the matrix. For `as.symmetric.matrix`, this is an
  optional argument that specifies the dimension of the symmetric
  matrix. This applies only when `x` is a vector. It defaults to `NULL`,
  in which case it calculates `k <- (-1 + sqrt(1 + 8 * length(x)))/ 2`.

- parm:

  This is a vector of parameters passed to the model specification.

- Data:

  This is the list of data passed to the model specification. For more
  information, see
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

- a, b:

  These optional arguments allow the elements of `x` to be bound to the
  interval \\\[a,b\]\\. For example, elements of a correlation matrix
  are in the interval \\\[-1,1\]\\.

- restrict:

  Logical. If `TRUE`, then `x[1,1]` is restricted to 1. This is useful
  in multinomial probit, for example. The variable, `luciferMatrix`, is
  created in a new environment, `LDEnv` so `as.parm.matrix` can keep
  track of changes from iteration to iteration.

- chol:

  Logical. If `TRUE`, then x is an upper-triangular matrix.

- Sigma:

  This is a covariance matrix, \\\Sigma\\, and may be entered either as
  a matrix or vector.

- Model:

  This is a model specification function. For more information, see
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

- Method:

  This accepts a quoted string. For `Hessian`, it defaults to
  `Method="Richardson"`, which uses Richardson extrapolation. For
  `Jacobian`, it defaults to `Method="simple"`, which uses
  finite-differencing. Richardson extrapolation is more accurate, but
  slower to calculate. Since error due to finite-differencing propagates
  through to higher derivatives, finite-differencing should not be used
  when approximating a Hessian matrix. Another method called automatic
  differentiation is not currently available here, but should be more
  accurate, though even slower to calculate. Another popular alternative
  is to use the
  [`BayesianBootstrap`](https://robustecologies.github.io/lucifer/reference/BayesianBootstrap.md)
  on the data. For `CovEstim`, this accepts `Method="Hessian"`,
  `Method="Identity"` (which simply assigns an identity matrix),
  `Method="OPG"` (which calculates the sum of outer products of
  record-level gradients), or `Method="Sandwich"`, which is the sandwich
  estimator and combines the Hessian and OPG estimates.

- N:

  This required argument accepts a positive integer that indicates the
  number of nodes.

- dims:

  This required argument indicates the dimension of the integral and
  accepts a positive integer.

- rule:

  This is an optional argument that accepts a univariate Gauss-Hermite
  quadrature rule. Usually, this argument is left empty. A rule may be
  supplied that differs from the traditional rule, such as when
  constraints have been observed, and one or more nodes and weights were
  adjusted.

- Interval:

  This accepts a small scalar number for precision.

- diag:

  Logical. If `TRUE`, then the elements in the main diagonal are also
  returned.

- file:

  This is the name of the file from which the numeric data matrix will
  be imported or read.

- header:

  Logical. When `TRUE`, the first row of the file must contain names of
  the columns, and will be converted to the column names of the numeric
  matrix. When `FALSE`, the first row of the file contains data, not
  column names.

- sep:

  This argument indicates a character with which it will separate fields
  when creating column vectors. For example, to read a comma-separated
  file (.csv), use `sep=","`.

- nrow:

  This is the number of rows of the numeric matrix, and defaults to
  `nrow=0`. If the number is known, the function will perform noticeably
  faster when it does not have to check.

- samples:

  This is the number of samples to take from the numeric matrix. When
  `samples=0`, sampling is not performed and the entire matrix is
  returned.

- size:

  This is the batch size to be used only when reading a numeric matrix
  that is larger than the available computer memory (RAM), and only when
  `samples` is greater than zero. Sampling of a big data matrix is
  performed by first determining the records to keep, and then reading
  batches, one by one, and keeping the matching records.

- na.rm:

  Logical. When `TRUE`, rows with missing values are removed from the
  matrix after it is read. Rather than removing missing values, the user
  may consider imputing missing values inside the model, or before the
  model with the
  [`MISS`](https://robustecologies.github.io/lucifer/reference/MISS.md)
  function. Examples of within-model imputation may be found in the
  accompanying "Examples" vignette.

- J:

  This required argument indicates the dimension of the integral and
  accepts a positive integer.

- K:

  This required argument indicates the accuracy and accepts a positive
  integer. Larger values result in many more integration nodes.

- theta.y:

  This accepts a vector of posterior samples of a discrete Markov chain,
  and defaults to `NULL`. If used, the order of the samples affects the
  transition probability.

- y.theta:

  This accepts a vector of data that are samples of a discrete
  distribution, and defaults to `NULL`. If used, the order of the
  samples affects the transition probability.

- p.theta:

  This accepts a matrix of prior probabilities for a transition matrix,
  and defaults to `NULL`. If used, each row must sum to 1.

## Value

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

## Details

Implementation of `as.indicator.matrix`. Refer to the package vignettes
and the cited references for a complete algorithmic and mathematical
description.

Implementation of `as.inverse`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `as.parm.matrix`. Refer to the package vignettes and
the cited references for a complete algorithmic and mathematical
description.

Implementation of `as.positive.definite`. Refer to the package vignettes
and the cited references for a complete algorithmic and mathematical
description.

Implementation of `as.positive.semidefinite`. Refer to the package
vignettes and the cited references for a complete algorithmic and
mathematical description.

Implementation of `as.symmetric.matrix`. Refer to the package vignettes
and the cited references for a complete algorithmic and mathematical
description.

Implementation of `Cov2Cor`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `CovEstim`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `GaussHermiteCubeRule`. Refer to the package vignettes
and the cited references for a complete algorithmic and mathematical
description.

Implementation of `Hessian`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `is.positive.definite`. Refer to the package vignettes
and the cited references for a complete algorithmic and mathematical
description.

Implementation of `is.positive.semidefinite`. Refer to the package
vignettes and the cited references for a complete algorithmic and
mathematical description.

Implementation of `is.square.matrix`. Refer to the package vignettes and
the cited references for a complete algorithmic and mathematical
description.

Implementation of `is.symmetric.matrix`. Refer to the package vignettes
and the cited references for a complete algorithmic and mathematical
description.

Implementation of `Jacobian`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `logdet`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `lower.triangle`. Refer to the package vignettes and
the cited references for a complete algorithmic and mathematical
description.

Implementation of `read.matrix`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `SparseGrid`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `TransitionMatrix`. Refer to the package vignettes and
the cited references for a complete algorithmic and mathematical
description.

Implementation of `tr`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `upper.triangle`. Refer to the package vignettes and
the cited references for a complete algorithmic and mathematical
description.

## References

Higham, N.J. (2002). "Computing the Nearest Correlation Matrix - a
Problem from Finance". *IMA Journal of Numerical Analysis*, 22, p.
329–343.
[doi:10.1093/imanum/22.3.329](https://doi.org/10.1093/imanum/22.3.329)

## See also

[`BayesianBootstrap`](https://robustecologies.github.io/lucifer/reference/BayesianBootstrap.md),
[`Cov2Prec`](https://robustecologies.github.io/lucifer/reference/Precision.md),
[`cov2cor`](https://rdrr.io/r/stats/cor.html),
[`ddirichlet`](https://robustecologies.github.io/lucifer/reference/dist.Dirichlet.md),
[`GaussHermiteQuadRule`](https://robustecologies.github.io/lucifer/reference/Math.md),
[`isSymmetric`](https://rdrr.io/r/base/isSymmetric.html),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`lower.tri`](https://rdrr.io/r/base/lower.tri.html),
[`MISS`](https://robustecologies.github.io/lucifer/reference/MISS.md),
[`Prec2Cov`](https://robustecologies.github.io/lucifer/reference/Precision.md),
[`solve`](https://rdrr.io/r/base/solve.html), and
[`upper.tri`](https://rdrr.io/r/base/lower.tri.html).

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate positive-definite matrix
A <- as.positive.definite(matrix(c(2, 1, 1, 1), 2, 2))

# Check positive semi-definiteness
is.positive.semidefinite(A)

# Create diagonal matrix
D <- as.parm.matrix(fit, a = 1, p = 3, restrict = TRUE)
} # }

if (FALSE) { # \dontrun{
## see package vignettes for full examples involving as.inverse
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving as.parm.matrix
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving as.positive.definite
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving as.positive.semidefinite
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving as.symmetric.matrix
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving Cov2Cor
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving CovEstim
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving GaussHermiteCubeRule
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving Hessian
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving is.positive.definite
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving is.positive.semidefinite
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving is.square.matrix
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving is.symmetric.matrix
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving Jacobian
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving logdet
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving lower.triangle
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving read.matrix
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving SparseGrid
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving TransitionMatrix
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving tr
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving upper.triangle
} # }
```
