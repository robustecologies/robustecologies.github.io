# Bayes' theorem

Bayes' theorem shows the relation between two conditional probabilities
that are the reverse of each other. This theorem is named after Reverend
Thomas Bayes (1702-1761), and is also referred to as Bayes' law or
Bayes' rule (Bayes and Price, 1763). Bayes' theorem expresses the
conditional probability, or "posterior probability", of an event \\A\\
after \\B\\ is observed in terms of the "prior probability" of \\A\\,
prior probability of \\B\\, and the conditional probability of \\B\\
given \\A\\. Bayes' theorem is valid in all common interpretations of
probability. This function provides one of several forms of calculations
that are possible with Bayes' theorem.

## Usage

``` r
BayesTheorem(PrA, PrBA)
```

## Arguments

- PrA:

  This required argument is the prior probability of \\A\\, or
  \\\Pr(A)\\.

- PrBA:

  This required argument is the conditional probability of \\B\\ given
  \\A\\ or \\\Pr(B \| A)\\, and is known as the data, evidence, or
  likelihood.

## Value

The `BayesTheorem` function returns the conditional probability of \\A\\
given \\B\\, known in Bayesian inference as the posterior. The returned
object is of class `bayestheorem`.

## Details

Bayes' theorem provides an expression for the conditional probability of
\\A\\ given \\B\\, which is equal to

\$\$\Pr(A \| B) = \frac{\Pr(B \| A)\Pr(A)}{\Pr(B)}\$\$

For example, suppose one asks the question: what is the probability of
going to Hell, conditional on consorting (or given that a person
consorts) with Lucifer. By replacing \\A\\ with \\Hell\\ and \\B\\ with
\\Consort\\, the question becomes

\$\$\Pr(\mathrm{Hell} \| \mathrm{Consort}) = \frac{\Pr(\mathrm{Consort}
\| \mathrm{Hell})\Pr(\mathrm{Hell})}{\Pr(\mathrm{Consort})}\$\$

Note that a common fallacy is to assume that \\\Pr(A \| B) = \Pr(B \|
A)\\, which is called the conditional probability fallacy.

Another way to state Bayes' theorem (and this is the form in the
provided function) is

\$\$\Pr(A_i \| B) = \frac{\Pr(B \| A_i)\Pr(A_i)}{\Pr(B \| A_i)\Pr(A_i)
+\dots+ \Pr(B \| A_n)\Pr(A_n)}\$\$

## References

Bayes, T. and Price, R. (1763). "An Essay Towards Solving a Problem in
the Doctrine of Chances". By the late Rev. Mr. Bayes, communicated by
Mr. Price, in a letter to John Canton, M.A. and F.R.S. *Philosophical
Transactions of the Royal Statistical Society of London*, 53, p.
370–418.

## See also

[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md), and
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Pr(Hell|Consort) =
PrA <- c(0.75,0.25)
PrBA <- c(6/9, 5/7)
BayesTheorem(PrA, PrBA)
} # }
```
