# Declarative probabilistic model specification

Compiles a probabilistic notation string into a lucifer-compatible Model
function, data template builder, and initial value generator. Supports
both plain text and LaTeX notation, reusing the parsing infrastructure
from
[`plot_dag`](https://robustecologies.github.io/lucifer/reference/plot_dag.md).

## Usage

``` r
model_spec(
  formula,
  data = NULL,
  monitor = "LP",
  compile = TRUE,
  backend = "auto"
)
```

## Arguments

- formula:

  Character string with one statement per line. Each line is either a
  stochastic statement (`LHS ~ Distribution(args)`), a deterministic
  assignment (`LHS = expression`), or a comment/blank line (ignored).
  LaTeX markup is stripped automatically, so notation copied from papers
  works directly.

- data:

  Optional character vector of variable names to force-classify as data,
  overriding automatic detection.

- monitor:

  Character vector of quantities to monitor (default `"LP"`).

- compile:

  Logical. If `FALSE`, returns the intermediate representation without
  compiling to a function.

- backend:

  Character string specifying the evaluation backend. One of `"r"` (pure
  R model function), `"cpp"` (compile to C++ via
  [`compile`](https://robustecologies.github.io/lucifer/reference/compile.md)),
  or `"auto"` (default, use C++ if Rcpp is available, otherwise fall
  back to R).

## Value

An S3 object of class `model_spec` containing:

- `$Model`:

  Compiled function(parm, Data) ready for
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

- `$data_template`:

  Function(...) that builds the Data list from user-provided data.

- `$initial_values`:

  Function(Data) that generates reasonable starting values.

- `$ir`:

  The intermediate representation (class `model_ir`).

- `$code`:

  Generated R source code as a character string.

- `$formula`:

  The original input text.

- `$parameters`:

  Data frame of parameter names, distributions, and constraints.

- `$data_requirements`:

  Data frame of required data fields.

## Details

**Two-phase compilation.** The system first parses the notation into a
structured intermediate representation (IR) that classifies every
variable as data, parameter, or hyperparameter, infers constraints from
distribution support, and validates consistency. The IR is then compiled
into R code following the exact patterns used in lucifer Model
functions: parameter extraction from the `parm` vector, constraint
application via
[`interval`](https://robustecologies.github.io/lucifer/reference/interval.md),
log-prior computation, deterministic relationships, log-likelihood,
posterior assembly, and the five-element return list.

**Parameterization convention.** Each distribution name maps to exactly
one parameterization to eliminate ambiguity. `Normal(mu, sigma)` uses
standard deviation; `NormalV(mu, var)` uses variance;
`NormalP(mu, prec)` uses precision. The full mapping table is documented
in the package vignette.

**Scope.** The DSL handles the majority of common Bayesian models: GLMs,
hierarchical models, multivariate models with standard priors and link
functions. Models requiring recursive time indexing, mixture allocation,
missing data augmentation, or complex control flow should be written
manually. The
[`code()`](https://robustecologies.github.io/lucifer/reference/code.md)
method provides a scaffold that can be copied and extended for such
cases.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`plot_dag`](https://robustecologies.github.io/lucifer/reference/plot_dag.md),
[`as.parm.names`](https://robustecologies.github.io/lucifer/reference/as.parm.names.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Linear regression
spec <- model_spec("
  y ~ Normal(mu, sigma)
  mu = X %*% beta
  beta[j] ~ Normal(0, 100), j = 1,...,J
  sigma ~ HalfCauchy(25)
")

# Inspect
print(spec)
code(spec)
plot(spec)

# Use with lucifer
N <- 100; J <- 3
X <- cbind(1, matrix(rnorm(N * (J - 1)), N, J - 1))
y <- X %*% c(2, -1, 0.5) + rnorm(N)
Data <- spec$data_template(y = y, X = X, J = J)
Initial.Values <- spec$initial_values(Data)
fit <- lucifer(spec$Model, Data, Initial.Values,
               Iterations = 10000, Algorithm = "NUTS")

# LaTeX notation (copy-pasted from paper)
spec2 <- model_spec("
  $\\textbf{y} \\sim \\mathcal{N}(\\mu, \\sigma)$
  $\\mu = \\textbf{X}\\beta$
  $\\beta_j \\sim \\mathcal{N}(0, 100), \\quad j=1,\\dots,J$
  $\\sigma \\sim \\mathcal{HC}(25)$
")
print(spec2)
} # }
```
