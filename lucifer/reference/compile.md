# Compile a model_spec to C++ for fast evaluation

Generates C++ source code from the model_spec intermediate
representation, compiles it via
[`Rcpp::sourceCpp()`](https://rdrr.io/pkg/Rcpp/man/sourceCpp.html), and
returns a compiled model object that can be passed directly to
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

## Usage

``` r
compile(spec, Data, cache = TRUE, verbose = TRUE)
```

## Arguments

- spec:

  A model_spec object.

- Data:

  The Data list (from `spec$data_template(...)`).

- cache:

  Logical; cache compiled objects for reuse. Default TRUE.

- verbose:

  Logical; print compilation messages. Default TRUE.

## Value

A `compiled_model` S3 object, or the original spec on failure.

## Details

The compiled model evaluates the log-posterior and gradient entirely in
C++ without R callbacks. For gradient-based samplers (NUTS, HMC, MALA),
this eliminates the 2n+1 R function calls per gradient evaluation and
provides 10-100x speedup.

The compilation uses
[`Rcpp::sourceCpp()`](https://rdrr.io/pkg/Rcpp/man/sourceCpp.html) and
requires a C++ compiler. Compiled objects are cached by model hash in
`tools::R_user_dir("lucifer", "cache")`.

If compilation fails, the function falls back to the R model with a
warning.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
spec <- model_spec("
    y ~ Normal(mu, sigma)
    mu = X %*% beta
    beta ~ Normal(0, 10)
    sigma ~ HalfCauchy(2.5)
")
Data <- spec$data_template(y = y, X = X)
compiled <- compile(spec, Data)
fit <- lucifer(compiled, Data, spec$initial_values(Data),
               Algorithm = "NUTS", Iterations = 5000)
} # }
```
