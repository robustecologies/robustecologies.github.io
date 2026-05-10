# Summary method for demonoid objects

The `summary.demonoid` function provides a structured summary of MCMC
results from
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
including acceptance rates, convergence assessment, DIC, effective
sample sizes, and posterior summary tables. For multi-chain runs, it
additionally reports per-chain acceptance rates, the Gelman-Rubin
diagnostic, split R-hat, and an overall convergence assessment.
Divergence counts from gradient-based samplers (HMC, NUTS, etc.) are
reported when present.

## Usage

``` r
# S3 method for class 'demonoid'
summary(object, ...)
```

## Arguments

- object:

  An object of class `demonoid`, as returned by
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

- ...:

  Additional arguments (currently unused).

## Value

Returns the input object invisibly. The summary is printed as a
side-effect.

## Details

Produces a tabular summary of an MCMC fit produced by
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).
Summary of the content is given below. The returned object has class
`summary.<class>` and carries marginal posterior quantiles, effective
sample sizes where applicable, and diagnostic flags. Printing the
summary object yields the human-readable table; subscripting it with
`` `$` `` exposes the underlying numeric matrix.

## References

Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., &
Rubin, D. B. (2013). *Bayesian Data Analysis* (3rd ed.). Chapman &
Hall/CRC. ISBN 9781439840955.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`print.demonoid`](https://robustecologies.github.io/lucifer/reference/print.demonoid.md),
[`Gelman.Diagnostic`](https://robustecologies.github.io/lucifer/reference/Gelman.Diagnostic.md),
[`Rhat`](https://robustecologies.github.io/lucifer/reference/Rhat.md),
[`ESS.bulk`](https://robustecologies.github.io/lucifer/reference/ESS.bulk.md),
[`ESS.tail`](https://robustecologies.github.io/lucifer/reference/ESS.tail.md),
[`Consort`](https://robustecologies.github.io/lucifer/reference/Consort.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## Single-chain example
fit <- lucifer(Model, MyData, Initial.Values,
    Iterations=10000, Algorithm="NUTS",
    Specs=list(A=500, delta=0.65, epsilon=0.01, Lmax=20))
summary(fit)

## Multi-chain example
fit <- lucifer(Model, MyData, Initial.Values,
    Iterations=10000, Algorithm="NUTS",
    Specs=list(A=500, delta=0.65, epsilon=0.01, Lmax=20),
    Chains=4, CPUs=4)
summary(fit)
} # }
```
