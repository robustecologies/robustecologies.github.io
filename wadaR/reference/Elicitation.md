# Prior elicitation

Prior elicitation is the act of inducing personal opinion to be
expressed by the probabilities the person associates with an event
(Savage, 1971). The `elicit` function elicits personal opinion and the
`delicit` function estimates probability density to be used with model
specification in the
[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md), or
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md)
functions.

See Details.

## Usage

``` r
delicit(theta, x, a = -Inf, b = Inf, log = FALSE)

elicit(n, cats, cat.names, show.plot = FALSE)
```

## Arguments

- theta:

  This is a scalar or vector of parameters for which the density is
  estimated with respect to the kernel density estimate of `x`.

- x:

  This is the elicited vector.

- a:

  This is an optional lower bound for support.

- b:

  This is an optional upper bound for support.

- log:

  Logical. If `log=TRUE`, then the logarithm of the density is returned.

- n:

  This is the number of chips.

- cats:

  This is a vector of \\k\\ categories, bins, or intervals. When the
  variable is continuous, the mid-point of each category is used. For
  example, if the continuous interval \[0,1\] has 5 equal-sized
  categories, then `cats=c(0.1,0.3,0.5,0.7,0.9)`.

- cat.names:

  This is a vector of category names. For example, if the continuous
  interval \[0,1\] has 5 equal-sized categories, then one way of naming
  the categories may be
  `cat.names=c("0:<.2", ".2:<.4", ".4:<.6", ".6:<.8", ".8:1")`.

- show.plot:

  Logical. If `show.plot=TRUE`, then a barplot is shown after each
  allocation of chips.

## Value

`delicit` returns a vector of density values (or log-density values if
`log=TRUE`).

`elicit` returns a vector of elicited values corresponding to the
allocated chips across categories.

See Details.

## Details

The `elicit` function elicits a univariate, discrete, non-conjugate,
informative, prior probability distribution by offering a number of
chips (specified as `n` by the statistician) for the user to allocate
into categories specified by the statistician. The results of multiple
elicitations (meaning, with multiple people), each the output of
`elicit`, may be combined with the `c` function in base R.

This discrete distribution is included with the data for a model and
supplied to a model specification function, where in turn it is supplied
to the `delicit` function, which estimates the density at the current
value of the prior distribution, \\p(\theta)\\. The prior distribution
may be either continuous or discrete, will be proper, and may have
bounded support (constrained to an interval).

For a minimal example, a statistician elicits the prior probability
distribution for a regression effect, \\\beta\\. Non-statisticians would
not be asked about expected parameters, but could be asked about how
much \\\textbf{y}\\ would be expected to change given a one-unit change
in \\\textbf{x}\\. After consulting with others who have prior
knowledge, the support does not need to be bounded, and their guesses at
the range result in the statistician creating 5 categories from the
interval \[-1,4\], where each interval has a width of one. The
statistician schedules time with 3 people, and each person participates
when the statistician runs the following R code:

`x <- elicit(n=10, cats=c(-0.5, 0.5, 1.5, 2.5, 3.5), cat.names=c("-1:<0", "0:<1", "1:<2", "2:<3", "3:4"), show.plot=TRUE)`

Each of the 3 participants receives 10 chips to allocate among the 5
categories according to personal beliefs in the probability of the
regression effect. When the statistician and each participant accept
their elicited distribution, all 3 vectors are combined into one vector.
In the model form, the prior is expressed as

\$\$p(\beta) \sim \mathcal{EL}\$\$

and the code for the model specification is

`elicit.prior <- delicit(beta, x, log=TRUE)`

This method is easily extended to priors that are multivariate,
correlated, or conditional.

As an alternative, Hahn (2006) also used a categorical approach,
eliciting judgements about the relative likelihood of each category, and
then minimizes the KLD (for more information on KLD, see the
[`KLD`](https://robustecologies.github.io/lucifer/reference/KLD.md)
function).

Implementation of `elicit`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

## References

Hahn, E.D. (2006). "Re-examining Informative Prior Elicitation Through
the Lens of Markov chain Monte Carlo Methods". *Journal of the Royal
Statistical Society*, A 169 (1), p. 37–48.
[doi:10.1111/j.1467-985X.2005.00381.x](https://doi.org/10.1111/j.1467-985X.2005.00381.x)

Savage, L.J. (1971). "Elicitation of Personal Probabilities and
Expectations". *Journal of the American Statistical Association*,
66(336), p. 783–801.
[doi:10.1080/01621459.1971.10482346](https://doi.org/10.1080/01621459.1971.10482346)

## See also

[`de.Finetti.Game`](https://robustecologies.github.io/lucifer/reference/de.Finetti.Game.md),
[`KLD`](https://robustecologies.github.io/lucifer/reference/KLD.md),
[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md), and
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
x <- c(1,2,2,3,3,3,4,7,8,8,9,10) #Elicited with elicit function
theta <- seq(from=-5,to=15,by=.1)
plot(theta, delicit(theta,x), type="l", xlab=expression(theta),
     ylab=expression("p(" * theta * ")"))
} # }

if (FALSE) { # \dontrun{
## see package vignettes for full examples involving elicit
} # }
```
