# Plot Hellinger distances

This function plots Hellinger distances in an object of class `bmk`.

## Usage

``` r
# S3 method for class 'bmk'
plot(
  x,
  col = colorRampPalette(c("black", "red"))(100),
  title = "",
  PDF = FALSE,
  Parms = NULL,
  ...
)
```

## Arguments

- x:

  This required argument is an object of class `bmk`. See the
  [`BMK.Diagnostic`](https://robustecologies.github.io/lucifer/reference/BMK.Diagnostic.md)
  function for more information.

- col:

  This argument specifies the colors of the cells. By default, the
  `colorRampPalette` function colors large Hellinger distances as `red`,
  small as `black`, and provides 100 color gradations.

- title:

  This argument specifies the title of the plot, and the default does
  not include a title.

- PDF:

  Logical. When `TRUE`, the plot is saved as a .pdf file.

- Parms:

  This argument accepts a vector of quoted strings to be matched for
  selecting parameters for plotting. This argument defaults to `NULL`
  and selects every parameter for plotting. Each quoted string is
  matched to one or more parameter names with the `grep` function. For
  example, if the user specifies `Parms=c("eta", "tau")`, and if the
  parameter names are beta\[1\], beta\[2\], eta\[1\], eta\[2\], and tau,
  then all parameters will be selected, because the string `eta` is
  within `beta`. Since `grep` is used, string matching uses regular
  expressions, so beware of meta-characters, though these are
  acceptable: ".", "\[", and "\]".

- ...:

  Additional arguments are unused.

## Value

See Details.

## Details

The `plot.bmk` function plots the Hellinger distances in an object of
class `bmk`. This is useful for quickly finding portions of chains with
large Hellinger distances, which indicates non-stationarity and
non-convergence.

## See also

[`BMK.Diagnostic`](https://robustecologies.github.io/lucifer/reference/BMK.Diagnostic.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(lucifer)
N <- 1000
J <- 10
Theta <- matrix(runif(N*J),N,J)
colnames(Theta) <- paste("beta[", 1:J, "]", sep="")
for (i in 2:N) {Theta[i,1] <- Theta[i-1,1] + rnorm(1)}
HD <- BMK.Diagnostic(Theta, batches=10)
plot(HD, title="Hellinger distance between batches")
} # }
```
