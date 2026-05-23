# Plot a numerical matrix

This function plots a numerical matrix, and is often used to plot the
following matrices: correlation, covariance, distance, and precision.

## Usage

``` r
plotMatrix(
  x,
  col = colorRampPalette(c("red", "black", "green"))(100),
  cex = 1,
  circle = TRUE,
  order = FALSE,
  zlim = NULL,
  title = "",
  Style = NULL,
  PDF = FALSE,
  ...
)
```

## Arguments

- x:

  A numerical matrix, or an object of class `bayesfactor`, `demonoid`,
  `iterquad`, `laplace`, `pmc`, `posteriorchecks`, or `vb`. See the
  details section for more information regarding these classes. One
  component of a blocked proposal covariance matrix must be pointed to
  explicitly, rather than to the object of class `demonoid`.

- col:

  The colors of the circles. By default, the `colorRampPalette` function
  colors strong positive correlation as `green`, zero correlation as
  `black`, and strong negative correlation as `red`, and provides 100
  color gradations.

- cex:

  When `circle=TRUE`, this specifies the size of the marginal text, the
  names of the parameters or variables, and defaults to 1.

- circle:

  Logical. When `TRUE`, each element in the numeric matrix is
  represented with a circle, and a larger circle is assigned to elements
  that are farther from zero. Also, when `TRUE`, the gradation scale
  does not appear to the right of the plot.

- order:

  Logical. Defaults to `FALSE`, and presents the parameters or variables
  in the same order as in the numeric matrix. When `TRUE`, the
  parameters or variables are ordered using principal components
  analysis.

- zlim:

  When `circle=FALSE`, the gradation scale may be constrained to an
  interval by `zlim`, such as `zlim=c(-1,1)`, and only values within the
  interval are plotted.

- title:

  The title of the plot. The default does not include a title. When `x`
  is of class `posteriorchecks`, the title is changed to
  `Posterior Correlation`.

- Style:

  Optional character string selecting the plot style: `"Circle"`
  (default) or `"Heatmap"`. When supplied, this overrides the `circle`
  argument.

- PDF:

  Logical. When `TRUE`, the plot is saved as a .pdf file.

- ...:

  Additional arguments are unused.

## Value

See Details.

## Details

The `plotMatrix` function produces one of two styles of plots, depending
on the `circle` argument. A \\K \times K\\ numeric matrix of \\K\\
parameters or variables is plotted. The plot is a matrix of the same
dimensions, in which each element is colored (and sized, when
`circle=TRUE`) according to its value.

Although `plotMatrix` does not provide the same detail as a numeric
matrix, it is easier to discover elements of interest according to color
(and size when `circle=TRUE`).

The `plotMatrix` function is not inherently Bayesian, and does not
include uncertainty in matrices. Nonetheless, it is included because it
is a useful graphical presentation of numeric matrices, and is
recommended to be used with the posterior correlation matrix in an
object of class `posteriorchecks`.

When `x` is an object of class `bayesfactor`, matrix `B` is plotted.
When `x` is an object of class `demonoid` (if it is a matrix),
`iterquad`, `laplace`, `pmc`, or `vb`, the covariance matrix `Covar` is
plotted. When `x` is an object of class `posteriorchecks`, the posterior
correlation matrix is plotted.

This is a modified version of the `circle.corr` function of Taiyun Wei.

## See also

[`PosteriorChecks`](https://robustecologies.github.io/lucifer/reference/PosteriorChecks.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(mtcars)
plotMatrix(cor(mtcars),
  col=colorRampPalette(c("green","gray10","red"))(100),
  cex=1, circle=FALSE, order=TRUE)
plotMatrix(cor(mtcars),
  col=colorRampPalette(c("green","gray10","red"))(100),
  cex=1, circle=TRUE, order=TRUE)
} # }
```
