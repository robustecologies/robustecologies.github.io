# PMC RAM estimate

This function estimates the random-access memory (RAM) required to
update a given model and data with
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md).

*Warning:* Unwise use of this function may crash a computer, so please
read the details below.

## Usage

``` r
PMC.RAM(Model, Data, Iterations, Thinning, M, N)
```

## Arguments

- Model:

  A model specification function. For more information, see
  [`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md).

- Data:

  A list of data. For more information, see
  [`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md).

- Iterations:

  The number of iterations for which
  [`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md)
  would update. For more information, see
  [`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md).

- Thinning:

  The amount of thinning applied to the samples in
  [`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md).
  For more information, see
  [`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md).

- M:

  The number of mixture components in
  [`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md).

- N:

  The number of samples in
  [`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md).

## Value

A list with several components. Each component is an estimate in MB for
an object:

- alpha:

  Estimated size in MB of the matrix of mixture probabilities by
  iteration.

- Covar:

  Estimated size in MB of the covariance matrix or matrices.

- Data:

  Estimated size in MB of the list of data.

- Deviance:

  Estimated size in MB of the deviance vector.

- Initial.Values:

  Estimated size in MB of the matrix or vector of initial values.

- LH:

  Estimated size in MB of the \\N \times T \times M\\ array `LH`.

- LP:

  Estimated size in MB of the \\N \times T \times M\\ array `LP`.

- Model:

  Estimated size in MB of the model specification function.

- Monitor:

  Estimated size in MB of the \\N \times J\\ matrix `Monitor`.

- Posterior1:

  Estimated size in MB of the \\N \times J \times T \times M\\ array.

- Posterior2:

  Estimated size in MB of the \\N \times J\\ matrix.

- Summary:

  Estimated size in MB of the summary table.

- W:

  Estimated size in MB of the matrix of importance weights.

- Total:

  Estimated total size in MB of RAM required.

## Details

The `PMC.RAM` function uses the
[`object.size`](https://rdrr.io/r/utils/object.size.html) function to
estimate the size in MB of RAM required to update in
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md) for
a given model and data, and for a number of iterations and specified
thinning. When RAM is exceeded, the computer will crash. This function
can be useful when trying to estimate how many samples and iterations to
update a model without crashing the computer. However, when estimating
the required RAM, `PMC.RAM` actually creates several large objects, such
as `post` (see below). If too many iterations are given as an argument
to `PMC.RAM`, for example, then it will crash the computer while trying
to estimate the required RAM.

The best way to use this function is as follows. First, prepare the
model specification and list of data. Second, observe how much RAM the
computer is using at the moment, as well as the maximum available RAM.
The majority of the difference of these two is the amount of RAM the
computer may dedicate to updating the model. Next, use this function
with a small number of iterations. Note the estimated RAM. Increase the
number of iterations, and again note the RAM. Continue to increase the
number of iterations until, say, arbitrarily within 90% of the
above-mentioned difference in RAM.

The computer operating system uses RAM, as does any other software
running at the moment. R is currently using RAM, and other functions in
the `lucifer` package, and any other package that is currently
activated, are using RAM. There are numerous small objects that are not
included in the returned list, that use RAM. For example, perplexity is
a small vector, etc.

A potentially large object that is not included is a matrix used for
estimating
[`LML`](https://robustecologies.github.io/lucifer/reference/LML.md).

## See also

[`BigData`](https://robustecologies.github.io/lucifer/reference/BigData.md),
[`LML`](https://robustecologies.github.io/lucifer/reference/LML.md),
[`object.size`](https://rdrr.io/r/utils/object.size.html), and
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving PMC.RAM
} # }
```
