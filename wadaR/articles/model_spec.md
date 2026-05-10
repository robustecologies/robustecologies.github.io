# Probabilistic model specification with model_spec()

## Introduction

The **lucifer** package requires users to write Model functions
manually: extract parameters from a `parm` vector using positional
indices, compute log-priors, build the log-likelihood, assemble the
log-posterior, and return a five-element list. This workflow is powerful
and flexible, but the boilerplate is verbose, error-prone, and
repetitive across models that differ only in their probabilistic
specification. The
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
function addresses this by providing a declarative front-end where users
write probabilistic notation and the system compiles it into a valid
Model function, without replacing the manual approach for models that
need it.

The design is a two-phase compiler: text or LaTeX notation is parsed
into an intermediate representation (IR), which is then compiled into R
code that follows the exact patterns established across 108+ example
models. The generated function is indistinguishable from hand-written
code at runtime; there is no interpretation overhead. Users who outgrow
the DSL can call `code(spec)` to extract the generated source and extend
it manually, which means the DSL is both a productivity tool and a
learning scaffold.

## Quick start

A complete workflow with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
compresses the usual four-step setup (Data list, parameter names, Model
function, initial values) into a single declaration followed by data
binding.

``` r

library(lucifer)

# Simulate data
set.seed(42)
N <- 100
x <- rnorm(N)
y <- 3 + 1.5 * x + rnorm(N, sd = 0.8)
X <- cbind(1, x)
J <- ncol(X)

# Specify the model declaratively
spec <- model_spec("
  y ~ Normal(mu, sigma)
  mu = X %*% beta
  beta[j] ~ NormalV(0, 1000), j = 1,...,J
  sigma ~ HalfCauchy(25)
")

# Inspect
print(spec)
code(spec)

# Bind data and generate initial values
Data <- spec$data_template(y = y, X = X, J = J)
Initial.Values <- spec$initial_values(Data)

# Fit with NUTS
set.seed(123)
fit <- lucifer(spec$Model, Data, Initial.Values,
               Iterations = 1000, Algorithm = "NUTS")
print(fit)
```

That is the entire workflow. The rest of this document explains every
component in detail, covers the full distribution registry, demonstrates
progressively complex examples, and honestly delineates what the DSL can
and cannot handle.

## Function signature and return value

``` r

model_spec(formula, data = NULL, monitor = "LP", compile = TRUE)
```

The `formula` argument is a character string containing one statement
per line in either plain text or LaTeX notation (or a mixture of both).
Blank lines and lines starting with `#` are ignored. The `data` argument
is an optional character vector that force-classifies specified
variables as data, overriding auto-detection. The `monitor` argument
specifies quantities to monitor (default `"LP"`). When
`compile = FALSE`, the function returns only the intermediate
representation without generating executable code, which is useful for
debugging the parser.

The return value is an S3 object of class `model_spec` with the
following components:

``` r

spec$Model             # compiled function(parm, Data) for lucifer()
spec$data_template     # function(...) that builds the Data list
spec$initial_values    # function(Data) that generates starting values
spec$ir                # intermediate representation (class model_ir)
spec$code              # generated R source code as character string
spec$formula           # original input text
spec$parameters        # data.frame: name, distribution, constraint, is_vector
spec$data_requirements # data.frame: name, role (response or covariate)
spec$monitor           # character vector
```

## S3 methods

The [`print()`](https://rdrr.io/r/base/print.html) method produces a
concise one-screen summary showing the model notation, parameter table,
and data requirements. The
[`summary()`](https://rdrr.io/r/base/summary.html) method extends this
with the full variable classification, generated code, and constraint
details. The [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
method delegates to
[`plot_dag()`](https://robustecologies.github.io/lucifer/reference/plot_dag.md)
to visualize the model as a directed acyclic graph. The
[`code()`](https://robustecologies.github.io/lucifer/reference/code.md)
generic extracts and prints the generated R source code.

``` r

print(spec)
summary(spec)
code(spec)
```

``` r

plot(spec)
```

## Notation grammar

Each line in the formula string is one of three types.

**Stochastic statements** declare a random variable and its
distribution: `LHS ~ Distribution(arg1, arg2, ...)`, with an optional
trailing index range `, j = 1,...,J` or `, j = 1:J`. The left-hand side
can be a plain name (`sigma`), a bracket-indexed name (`beta[j]`), or an
underscore-indexed name (`beta_j`, restricted to single-letter
subscripts to avoid ambiguity with compound names like `mu_beta`).

**Deterministic statements** declare a computed quantity:
`LHS = R_expression`. The right-hand side is any valid R expression; it
is transcribed verbatim into the generated code, with data variables
automatically prefixed by `Data$`.

**Comments and blanks** (lines starting with `#` or empty) are ignored.

``` r

# All three line types
spec <- model_spec("
  # This is a comment
  y ~ Normal(mu, sigma)
  mu = X %*% beta
  beta[j] ~ NormalV(0, 1000), j = 1,...,J
  sigma ~ HalfCauchy(25)
")
```

## LaTeX notation

The parser strips LaTeX markup before processing, which means users can
paste notation directly from papers. The following transformations are
applied in order: dollar signs and `\textbf{}`, `\text{}`, `\mathcal{}`
wrappers are removed (content preserved); `\sim` becomes `~`;
superscripts `^{...}` are removed; Greek letter commands (`\alpha`,
`\beta`, `\sigma`, `\mu`, etc.) become their plain-text equivalents;
`\dots`, `\cdots`, `\ldots` become `...`; `\quad` and `\qquad` become
spaces; remaining backslash commands are removed; and stray braces
[`{}`](https://rdrr.io/r/base/Paren.html) are cleaned up.

``` r

spec_latex <- model_spec("
  $\\textbf{y} \\sim \\mathcal{N}(\\mu, \\sigma)$
  $\\mu = \\textbf{X}\\beta$
  $\\beta_j \\sim \\mathcal{N}_v(0, 1000), \\quad j=1,\\dots,J$
  $\\sigma \\sim \\mathcal{HC}(25)$
")

# Produces identical compiled code to the plain-text version
```

Note that R string escaping requires doubling backslashes: `\\sim` in an
R string becomes `\sim` in the notation. When reading from a file via
[`readLines()`](https://rdrr.io/r/base/readLines.html), this escaping is
not needed.

## Distribution registry

The DSL recognizes 86 distributions. Each has a canonical DSL name, zero
or more aliases (LaTeX abbreviations and lucifer function names), and a
fixed parameterization. Users can refer to any distribution by its
canonical name, any alias, or its lucifer density function name.

### Parameterization convention

The parameterization is explicit and unambiguous. There are no
heuristics: `Normal(mu, sigma)` always means standard deviation,
`NormalV(mu, var)` always means variance, and `NormalP(mu, prec)` always
means precision. This design sacrifices LaTeX purity for correctness,
since silent parameterization errors in MCMC would be far worse than
requiring users to specify `NormalV` when they mean variance.

### Naming styles

Every distribution can be referenced in three ways:

``` r

# All three resolve to the same distribution
spec1 <- model_spec("y ~ Normal(mu, sigma)\nmu ~ Normal(0, 100)\nsigma ~ HalfCauchy(25)")
spec2 <- model_spec("y ~ N(mu, sigma)\nmu ~ N(0, 100)\nsigma ~ HC(25)")
spec3 <- model_spec("y ~ dnorm(mu, sigma)\nmu ~ dnorm(0, 100)\nsigma ~ dhalfcauchy(25)")
```

### Complete distribution table

The following table lists every distribution in the registry with its
canonical DSL name, lucifer density function, and parameterization.

| DSL name             | lucifer dfun    | Parameterization          | Constraint |
|:---------------------|:----------------|:--------------------------|:-----------|
| Normal               | dnorm           | mean, sd                  | none       |
| NormalV              | dnormv          | mean, var                 | none       |
| NormalP              | dnormp          | mean, prec                | none       |
| StudentT             | dst             | mu, sigma, nu             | none       |
| StudentTP            | dstp            | mu, tau, nu               | none       |
| Laplace              | dlaplace        | location, scale           | none       |
| LaplaceP             | dlaplacep       | mu, tau                   | none       |
| AsymmetricLaplace    | dalaplace       | location, scale, kappa    | none       |
| LogLaplace           | dllaplace       | location, scale           | positive   |
| PowerExponential     | dpe             | mu, sigma, kappa          | none       |
| NormalLaplace        | dnormlaplace    | mu, sigma, alpha, beta    | none       |
| ExGaussian           | dexgaussian     | mu, sigma, lambda         | none       |
| Stable               | dstable         | alpha, beta, gamma, delta | none       |
| VonMises             | dvonmises       | mu, kappa                 | none       |
| Simplex              | dsimplex        | mu, sigma                 | unit       |
| HalfCauchy           | dhalfcauchy     | scale                     | positive   |
| HalfNormal           | dhalfnorm       | scale                     | positive   |
| HalfT                | dhalft          | scale, nu                 | positive   |
| Gamma                | dgamma          | shape, rate               | positive   |
| InvGamma             | dinvgamma       | shape, scale              | positive   |
| Exponential          | dexp            | rate                      | positive   |
| Weibull              | dweibull        | shape, scale              | positive   |
| LogNormal            | dlnorm          | meanlog, sdlog            | positive   |
| LogNormalP           | dlnormp         | mu, tau                   | positive   |
| Pareto               | dpareto         | alpha                     | positive   |
| InvChiSq             | dinvchisq       | df, scale                 | positive   |
| InvGaussian          | dinvgaussian    | mu, lambda                | positive   |
| InvBeta              | dinvbeta        | a, b                      | positive   |
| GPD                  | dgpd            | mu, sigma, xi             | none       |
| HorseShoe            | dhs             | lambda, tau               | positive   |
| HuangWand            | dhuangwand      | nu, a, A                  | positive   |
| Lasso                | dlasso          | sigma, tau, lambda, a, b  | positive   |
| YangBerger           | dyangberger     | (none)                    | positive   |
| Tweedie              | dtweedie        | mu, phi, power            | positive   |
| HuangWandC           | dhuangwandc     | nu, a, A                  | positive   |
| YangBergerC          | dyangbergerc    | (none)                    | positive   |
| HyperG               | dhyperg         | g, alpha                  | positive   |
| Beta                 | dbeta           | shape1, shape2            | unit       |
| Uniform              | dunif           | min, max                  | bounded    |
| Poisson              | dpois           | lambda                    | none       |
| GenPoisson           | dgpois          | lambda, omega             | none       |
| Bernoulli            | dbern           | prob                      | none       |
| Binomial             | dbinom          | size, prob                | none       |
| NegBinomial          | dnbinom         | size, prob                | none       |
| Categorical          | dcat            | p                         | none       |
| MVN                  | dmvn            | mu, Sigma                 | none       |
| MVNC                 | dmvnc           | mu, U                     | none       |
| MVNP                 | dmvnp           | mu, Omega                 | none       |
| MVNPC                | dmvnpc          | mu, U                     | none       |
| MVT                  | dmvt            | mu, S, df                 | none       |
| MVTP                 | dmvtp           | mu, Omega, nu             | none       |
| MVTC                 | dmvtc           | mu, U, df                 | none       |
| MVTPC                | dmvtpc          | mu, U, nu                 | none       |
| MVCauchy             | dmvc            | mu, S                     | none       |
| MVCauchyC            | dmvcc           | mu, U                     | none       |
| MVCauchyP            | dmvcp           | mu, Omega                 | none       |
| MVCauchyPC           | dmvcpc          | mu, U                     | none       |
| MVLaplace            | dmvl            | mu, Sigma                 | none       |
| MVLaplaceC           | dmvlc           | mu, U                     | none       |
| MVPE                 | dmvpe           | mu, Sigma, kappa          | none       |
| MVPEC                | dmvpec          | mu, U, kappa              | none       |
| MVPolya              | dmvpolya        | alpha                     | none       |
| AML                  | daml            | mu, Sigma                 | none       |
| NormalMixture        | dnormm          | p, mu, sigma              | none       |
| LaplaceMixture       | dlaplacem       | p, location, scale        | none       |
| SkewLaplace          | dslaplace       | mu, alpha, beta           | none       |
| SkewDiscreteLaplace  | dsdlaplace      | p, q                      | none       |
| AsymmetricLogLaplace | dallaplace      | location, scale, kappa    | none       |
| Wishart              | dwishart        | nu, S                     | none       |
| WishartC             | dwishartc       | nu, S                     | none       |
| InvWishart           | dinvwishart     | nu, S                     | none       |
| InvWishartC          | dinvwishartc    | nu, S                     | none       |
| Dirichlet            | ddirichlet      | alpha                     | none       |
| CRMRF                | dcrmrf          | alpha, Omega              | none       |
| MatrixNormal         | dmatrixnorm     | M, U, V                   | none       |
| MatrixGamma          | dmatrixgamma    | alpha, beta, Sigma        | none       |
| InvMatrixGamma       | dinvmatrixgamma | alpha, beta, Psi          | none       |
| MatrixT              | dmatrixt        | M, Omega, Sigma, nu       | none       |
| NormalInvWishart     | dnorminvwishart | mu0, lambda, Sigma, S, nu | none       |
| NormalWishart        | dnormwishart    | mu0, lambda, Omega, S, nu | none       |
| SIW                  | dsiw            | nu, S, zeta, mu, delta    | none       |
| Zellner              | dzellner        | g, sigma, X               | none       |
| Truncated            | dtrunc          | spec, a, b                | none       |
| Clayton              | dclayton        | u, v, theta               | none       |
| Frank                | dfrank          | u, v, theta               | none       |
| CAR                  | dcar            | W, tau                    | none       |

Complete distribution registry {.table}

### Constraint inference

Constraints are inferred automatically from the distribution’s support:

- **Positive** (e.g., HalfCauchy, Gamma, Exponential, Weibull,
  LogNormal): the parameter is constrained via
  `interval(x, 1e-100, Inf)`.
- **Unit** (e.g., Beta, Simplex): the parameter is constrained via
  `interval(x, 1e-10, 1 - 1e-10)`.
- **Bounded** (Uniform): the parameter is constrained to the specified
  bounds.
- **None** (e.g., Normal, StudentT, Laplace): no constraint applied.

These constraints are applied in the generated code after parameter
extraction and before log-prior computation, exactly as in hand-written
Model functions.

## Variable classification

The DSL automatically classifies variables as data, parameters,
hyperpriors, or deterministic quantities using the following
priority-ordered rules:

1.  Variables explicitly listed in the `data` argument are classified as
    data.
2.  The left-hand side of a stochastic statement (`~`) is a parameter.
3.  Among stochastic LHS variables, the one whose value is not defined
    by a prior elsewhere (i.e., it appears as a LHS of `~` but no other
    stochastic statement provides its distribution) is the **likelihood
    response**; its LHS is data.
4.  Variables that appear only on the right-hand side and are not
    parameters are data.
5.  Parameters whose values appear as arguments in other parameters’
    priors are classified as hyperpriors.

In practice: given `y ~ Normal(mu, sigma)` and `mu ~ Normal(0, 100)` and
`sigma ~ HalfCauchy(25)`, the system recognizes that `mu` and `sigma`
have priors, `y` does not, so `y` is data. Any variable in deterministic
expressions that is not a parameter (like `X` in `mu = X %*% beta`) is
also data.

When auto-detection fails or is ambiguous, use the `data` argument to
override:

``` r

spec <- model_spec("
  y ~ Normal(mu, sigma)
  mu ~ Normal(0, 100)
  sigma ~ HalfCauchy(25)
", data = c("y"))
```

## Data template and initial values

The `data_template()` closure generates a properly structured Data list
for
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md).
It creates `parm.names` using
[`as.parm.names()`](https://robustecologies.github.io/lucifer/reference/as.parm.names.md),
adds position indices (`Data$pos.beta`, `Data$pos.sigma`, etc.), and
attaches the user-supplied data.

``` r

# Rebuild the linear regression spec for this demonstration
spec_reg <- model_spec("
  y ~ Normal(mu, sigma)
  mu = X %*% beta
  beta[j] ~ NormalV(0, 1000), j = 1,...,J
  sigma ~ HalfCauchy(25)
")

Data <- spec_reg$data_template(y = y, X = X, J = J)
str(Data)
```

The `initial_values()` closure generates constraint-aware starting
values: unconstrained parameters draw from `rnorm(n, 0, 0.1)`, positive
parameters from `runif(n, 0.1, 2)`, unit parameters from
`runif(n, 0.2, 0.8)`, and bounded parameters from
`runif(n, 0.3, 0.7) * (upper - lower) + lower`.

``` r

set.seed(1)
Initial.Values <- spec_reg$initial_values(Data)
print(Initial.Values)
```

## Examples

The following examples demonstrate progressively complex models,
comparing the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
declarative approach against the traditional manual Model function. In
every case, both approaches produce identical posterior inference when
run with the same random seed and algorithm.

### Example 1: simple normal model

The simplest case: estimate the mean and standard deviation of normally
distributed data.

**Traditional approach:**

``` r

set.seed(42)
y <- rnorm(50, mean = 5, sd = 2)

# Manual Model function
Model_manual <- function(parm, Data) {
    mu <- parm[Data$pos.mu]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    mu.prior <- dnorm(mu, 0, 100, log = TRUE)
    sigma.prior <- dhalfcauchy(sigma, 25, log = TRUE)
    LL <- sum(dnorm(Data$y, mu, sigma, log = TRUE))
    LP <- LL + mu.prior + sigma.prior
    yhat <- rnorm(length(Data$y), mu, sigma)
    list(LP = LP, Dev = -2 * LL, Monitor = LP,
         yhat = yhat, parm = parm)
}

# Manual Data list
parm.names <- as.parm.names(list(mu = 0, sigma = 0))
Data_manual <- list(
    y = y, N = length(y),
    parm.names = parm.names, mon.names = "LP",
    pos.mu = 1, pos.sigma = 2
)
Initial.Values <- c(0, 1)

set.seed(123)
fit_manual <- lucifer(Model_manual, Data_manual, Initial.Values,
                      Iterations = 1000, Algorithm = "NUTS")
```

**model_spec approach:**

``` r

spec <- model_spec("
  y ~ Normal(mu, sigma)
  mu ~ Normal(0, 100)
  sigma ~ HalfCauchy(25)
")

Data_spec <- spec$data_template(y = y)
Initial.Values <- spec$initial_values(Data_spec)

set.seed(123)
fit_spec <- lucifer(spec$Model, Data_spec, Initial.Values,
                    Iterations = 1000, Algorithm = "NUTS")
```

Both fits converge to the same posterior: \\\mu \approx 5\\, \\\sigma
\approx 2\\.

### Example 2: linear regression

``` r

set.seed(42)
N <- 200
X <- cbind(1, rnorm(N), rnorm(N))
J <- ncol(X)
beta_true <- c(2, -1.5, 0.8)
sigma_true <- 1.2
y <- X %*% beta_true + rnorm(N, 0, sigma_true)

spec <- model_spec("
  y ~ Normal(mu, sigma)
  mu = X %*% beta
  beta[j] ~ NormalV(0, 1000), j = 1,...,J
  sigma ~ HalfCauchy(25)
")

Data <- spec$data_template(y = y, X = X, J = J)
Initial.Values <- spec$initial_values(Data)

set.seed(123)
fit <- lucifer(spec$Model, Data, Initial.Values,
               Iterations = 1000, Algorithm = "NUTS")
```

The generated code (visible via `code(spec)`) follows the same pattern
as Examples.Rmd line 4357, using `X %*% beta` for the linear predictor.

### Example 3: logistic regression

``` r

set.seed(42)
N <- 300
X <- cbind(1, rnorm(N))
J <- ncol(X)
beta_true <- c(-0.5, 1.2)
p_true <- 1 / (1 + exp(-(X %*% beta_true)))
y <- rbinom(N, 1, p_true)

spec <- model_spec("
  y ~ Bernoulli(prob)
  prob = 1 / (1 + exp(-(X %*% beta)))
  beta[j] ~ Normal(0, 10), j = 1,...,J
")

Data <- spec$data_template(y = y, X = X, J = J)
Initial.Values <- spec$initial_values(Data)

set.seed(123)
fit <- lucifer(spec$Model, Data, Initial.Values,
               Iterations = 1000, Algorithm = "NUTS")
```

The DSL handles the inverse-logit link function as a standard
deterministic statement. Any valid R expression can appear on the
right-hand side of `=`.

### Example 4: Poisson regression

``` r

set.seed(42)
N <- 200
X <- cbind(1, rnorm(N))
J <- ncol(X)
beta_true <- c(1.5, 0.5)
lambda_true <- exp(X %*% beta_true)
y <- rpois(N, lambda_true)

spec <- model_spec("
  y ~ Poisson(lambda)
  lambda = exp(X %*% beta)
  beta[j] ~ Normal(0, 10), j = 1,...,J
")

Data <- spec$data_template(y = y, X = X, J = J)
Initial.Values <- spec$initial_values(Data)

set.seed(123)
fit <- lucifer(spec$Model, Data, Initial.Values,
               Iterations = 1000, Algorithm = "NUTS")
```

### Example 5: hierarchical linear model

This example demonstrates hyperprior specification, where the prior on
group-level parameters is itself governed by a hyperparameter.

``` r

set.seed(42)
N <- 200
G <- 5    # number of groups
group <- sample(1:G, N, replace = TRUE)
X <- cbind(1, rnorm(N))
J <- ncol(X)
tau_true <- 0.5
beta_true <- matrix(rnorm(G * J, 0, tau_true), G, J)
beta_true[, 1] <- beta_true[, 1] + 3
y <- rowSums(X * beta_true[group, ]) + rnorm(N, 0, 0.5)

spec <- model_spec("
  y ~ Normal(mu, sigma)
  mu = X %*% beta
  beta[j] ~ Normal(0, tau), j = 1,...,J
  tau ~ HalfCauchy(5)
  sigma ~ HalfCauchy(5)
")

Data <- spec$data_template(y = y, X = X, J = J)
Initial.Values <- spec$initial_values(Data)

set.seed(123)
fit <- lucifer(spec$Model, Data, Initial.Values,
               Iterations = 1000, Algorithm = "NUTS")
```

Here `tau` is classified as a hyperprior because it appears as an
argument in the prior for `beta`. The compiler sorts the code so that
`tau` is extracted and constrained before it is used in `beta.prior`.

### Example 6: beta regression

``` r

set.seed(42)
N <- 150
x <- rnorm(N)
mu_true <- 1 / (1 + exp(-(0.5 + 0.8 * x)))
phi_true <- 20
a_true <- mu_true * phi_true
b_true <- (1 - mu_true) * phi_true
y <- rbeta(N, a_true, b_true)

spec <- model_spec("
  y ~ Beta(a, b)
  a = mu * phi
  b = (1 - mu) * phi
  mu = 1 / (1 + exp(-(alpha + beta_x * x)))
  alpha ~ Normal(0, 10)
  beta_x ~ Normal(0, 10)
  phi ~ Gamma(0.1, 0.1)
")

Data <- spec$data_template(y = y, x = x)
Initial.Values <- spec$initial_values(Data)

set.seed(123)
fit <- lucifer(spec$Model, Data, Initial.Values,
               Iterations = 1000, Algorithm = "NUTS")
```

This demonstrates multiple deterministic statements and the ability to
express reparameterized likelihoods. The variables `a`, `b`, and `mu`
are all deterministic transformations.

### Example 7: using lucifer function names

Users familiar with lucifer’s native distribution functions can use them
directly:

``` r

spec <- model_spec("
  y ~ dnorm(mu, sigma)
  mu = X %*% beta
  beta[j] ~ dnormv(0, 1000), j = 1,...,J
  sigma ~ dhalfcauchy(25)
")
```

This is fully equivalent to using the canonical DSL names. Every density
function in lucifer (86 total) is registered as an alias and can be used
in the notation.

### Example 8: mixed notation styles

Plain text, LaTeX, and lucifer function names can be freely mixed within
the same specification:

``` r

spec <- model_spec("
  y ~ Normal(mu, sigma)
  mu = X %*% beta
  beta[j] ~ dnormv(0, 1000), j = 1,...,J
  sigma ~ dhalfcauchy(25)
")
```

### Example 9: passing model_spec directly to lucifer()

The
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
function recognizes `model_spec` objects and auto-unwraps them:

``` r

spec <- model_spec("
  y ~ Normal(mu, sigma)
  mu ~ Normal(0, 100)
  sigma ~ HalfCauchy(25)
")

Data <- spec$data_template(y = y)
Initial.Values <- spec$initial_values(Data)

# Pass spec directly instead of spec$Model
set.seed(123)
fit <- lucifer(spec, Data, Initial.Values,
               Iterations = 1000, Algorithm = "NUTS")
```

### Example 10: DAG visualization

The [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method
generates a directed acyclic graph showing the dependency structure:

``` r

spec <- model_spec("
  y ~ Normal(mu, sigma)
  mu = X %*% beta
  beta[j] ~ NormalV(0, 1000), j = 1,...,J
  sigma ~ HalfCauchy(25)
")

plot(spec)
```

This delegates to
[`plot_dag()`](https://robustecologies.github.io/lucifer/reference/plot_dag.md),
which also accepts `model_spec` objects directly:

``` r

plot_dag(spec)
```

## Code generation details

The
[`code()`](https://robustecologies.github.io/lucifer/reference/code.md)
method reveals the exact R code that the compiler produces.
Understanding this output is useful for debugging, learning lucifer’s
Model function conventions, and as a scaffold for extending models
beyond the DSL’s scope.

``` r

spec <- model_spec("
  y ~ Normal(mu, sigma)
  mu = X %*% beta
  beta[j] ~ NormalV(0, 1000), j = 1,...,J
  sigma ~ HalfCauchy(25)
")

code(spec)
```

The output follows this structure:

``` r

function(parm, Data) {
    # Parameter extraction
    beta <- parm[Data$pos.beta]
    sigma <- parm[Data$pos.sigma]

    # Constraints
    sigma <- interval(sigma, 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma

    # Log-priors
    beta.prior <- sum(dnormv(beta, 0, 1000, log = TRUE))
    sigma.prior <- dhalfcauchy(sigma, 25, log = TRUE)

    # Deterministic relationships
    mu <- Data$X %*% beta

    # Log-likelihood
    LL <- sum(dnorm(Data$y, mu, sigma, log = TRUE))

    # Log-posterior
    LP <- LL + beta.prior + sigma.prior

    # Posterior predictive
    yhat <- rnorm(length(Data$y), mu, sigma)

    list(LP = LP, Dev = -2 * LL, Monitor = LP,
         yhat = yhat, parm = parm)
}
```

The key conventions: parameters are extracted from `parm` via position
indices stored in `Data$pos.*`; constraints are applied immediately and
written back to `parm`; priors are computed as log-densities and summed
for vectorized parameters; data variables are accessed via `Data$`; the
log-posterior is the sum of log-likelihood and all log-priors; posterior
predictive samples (`yhat`) are drawn from the likelihood distribution;
and the return list contains `LP`, `Dev` (\\-2 \times \text{LL}\\),
`Monitor`, `yhat`, and `parm`.

## Compile-only mode

Setting `compile = FALSE` returns only the intermediate representation,
without generating executable code. This is useful for inspecting the
parser output:

``` r

ir_only <- model_spec("
  y ~ Normal(mu, sigma)
  mu ~ Normal(0, 100)
  sigma ~ HalfCauchy(25)
", compile = FALSE)

# No Model function available
is.null(ir_only$Model)
#> [1] TRUE

# But the IR is accessible
str(ir_only$ir)
```

## Scope and limitations

The DSL handles approximately 80% of models that users typically write
for lucifer. The remaining 20% require manual Model functions, and the
system is designed to be honest about this boundary.

### What model_spec handles well

- **Generalized linear models**: Normal, Poisson, Bernoulli, Binomial,
  Beta, Gamma, Negative Binomial regression with arbitrary link
  functions expressed as deterministic statements.
- **Hierarchical and multilevel models** with nested priors and
  hyperpriors.
- **Scalar, vector, and matrix parameters** with explicit indexing.
- **All 86 distributions** in the lucifer ecosystem, including
  multivariate (MVN, MVT, Wishart), mixture (NormalMixture), copula
  (Clayton, Frank), and spatial (CAR) distributions.
- **Standard linear predictors**: `mu = X %*% beta`,
  `mu = X %*% beta + Z %*% gamma`, etc.
- **Arbitrary deterministic transformations**: any valid R expression.

### What model_spec cannot handle

These patterns genuinely require a programming language and cannot be
expressed in a notation language. The DSL detects several of these and
issues warnings, suggesting manual Model writing.

**Time series with recursive structure:** AR, GARCH, and state-space
models require loops with lag indexing (e.g.,
`y[t] ~ Normal(phi * y[t-1], sigma)`). The DSL cannot express temporal
recursion.

**Mixture models with discrete latent allocation:** models like
`z[i] ~ Categorical(pi)` followed by
`y[i] ~ Normal(mu[z[i]], sigma[z[i]])` require array indexing by a
discrete latent variable.

**Missing data augmentation:** treating missing observations as
parameters that are sampled alongside other parameters requires
augmenting the `parm` vector with imputed values.

**Custom constraint patterns:** ordinal thresholds with
[`sort()`](https://rdrr.io/r/base/sort.html), sum-to-zero constraints,
or correlation matrices with Cholesky decomposition.

**Non-standard likelihoods:** CPO, ABC approximations, custom distance
functions, or likelihoods that require loops over observations.

**Complex control flow:** `if`/`else` statements within the model,
observation-level branching.

While the DSL can syntactically parse these patterns, the generated code
will be incorrect. For instance, a time series model compiles but
produces code referencing undefined index variables:

``` r

# This compiles but the generated code won't work correctly
spec_ts <- model_spec("
  y[t] ~ Normal(phi * y[t-1], sigma)
  phi ~ Normal(0, 1)
  sigma ~ HalfCauchy(5)
")
# The generated code shows Data$y[t-1] with undefined 't':
code(spec_ts)
```

The correct approach for such models is to write the Model function
manually with an explicit loop.

### The code() escape hatch

For models near the boundary of the DSL’s scope, the recommended
workflow is to use
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
to generate a scaffold, then extend the generated code manually:

``` r

# Step 1: get the scaffold
spec <- model_spec("
  y ~ Normal(mu, sigma)
  mu = X %*% beta
  beta[j] ~ NormalV(0, 1000), j = 1,...,J
  sigma ~ HalfCauchy(25)
")

# Step 2: extract the generated code
code(spec)

# Step 3: copy, paste, and modify for your specific needs
# For instance, add AR(1) structure:
Model_custom <- function(parm, Data) {
    beta <- parm[Data$pos.beta]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    phi <- parm[Data$pos.phi]  # Added manually

    beta.prior <- sum(dnormv(beta, 0, 1000, log = TRUE))
    sigma.prior <- dhalfcauchy(sigma, 25, log = TRUE)
    phi.prior <- dnorm(phi, 0, 1, log = TRUE)  # Added manually

    mu <- Data$X %*% beta
    # AR(1) residual structure added manually
    resid <- Data$y - mu
    LL <- dnorm(resid[1], 0, sigma, log = TRUE)
    for (t in 2:Data$N) {
        LL <- LL + dnorm(resid[t], phi * resid[t - 1], sigma, log = TRUE)
    }

    LP <- LL + beta.prior + sigma.prior + phi.prior
    yhat <- rnorm(Data$N, mu, sigma)
    list(LP = LP, Dev = -2 * LL, Monitor = LP,
         yhat = yhat, parm = parm)
}
```

This workflow preserves the DSL’s value as a starting point while
allowing full manual control.

## Comparison: traditional vs. model_spec

The following comparison uses identical random seeds to demonstrate that
both approaches produce the same posterior.

``` r

library(lucifer)
set.seed(42)

# Simulate data
N <- 100
x <- rnorm(N)
y <- 3 + 1.5 * x + rnorm(N, sd = 0.8)
X <- cbind(1, x)
J <- ncol(X)

# ========================================
# Traditional approach
# ========================================
Model_trad <- function(parm, Data) {
    beta <- parm[Data$pos.beta]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma

    beta.prior <- sum(dnormv(beta, 0, 1000, log = TRUE))
    sigma.prior <- dhalfcauchy(sigma, 25, log = TRUE)

    mu <- Data$X %*% beta
    LL <- sum(dnorm(Data$y, mu, sigma, log = TRUE))
    LP <- LL + beta.prior + sigma.prior

    yhat <- rnorm(length(Data$y), mu, sigma)
    list(LP = LP, Dev = -2 * LL, Monitor = LP,
         yhat = yhat, parm = parm)
}

parm.names <- as.parm.names(list(beta = rep(0, J), sigma = 0))
Data_trad <- list(
    y = y, X = X, J = J, N = N,
    parm.names = parm.names,
    mon.names = "LP",
    pos.beta = 1:J,
    pos.sigma = J + 1
)
IV_trad <- c(rep(0, J), 1)

set.seed(123)
fit_trad <- lucifer(Model_trad, Data_trad, IV_trad,
                    Iterations = 1000, Algorithm = "NUTS")

# ========================================
# model_spec approach
# ========================================
spec <- model_spec("
  y ~ Normal(mu, sigma)
  mu = X %*% beta
  beta[j] ~ NormalV(0, 1000), j = 1,...,J
  sigma ~ HalfCauchy(25)
")

Data_spec <- spec$data_template(y = y, X = X, J = J)
IV_spec <- spec$initial_values(Data_spec)

# Use same initial values for fair comparison
IV_spec <- IV_trad

set.seed(123)
fit_spec <- lucifer(spec$Model, Data_spec, IV_spec,
                    Iterations = 1000, Algorithm = "NUTS")

# Compare posteriors
summary_trad <- fit_trad$Summary1
summary_spec <- fit_spec$Summary1
cat("Traditional beta[1]:", summary_trad["beta[1]", "Mean"], "\n")
cat("model_spec beta[1]:", summary_spec["beta[1]", "Mean"], "\n")
# These will be identical (same seed, same algorithm, same code)
```

## Workflow recommendations

1.  **Start declarative, extend manually.** Use
    [`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
    for rapid prototyping. If the model compiles and runs, you are done.
    If it cannot express your model, call `code(spec)` and build from
    the scaffold.

2.  **Inspect the generated code.** Always call `code(spec)` at least
    once to verify the compilation matches your intent. The generated
    code is human-readable and follows established conventions.

3.  **Use `data_template()` for data management.** The template
    generates correct `parm.names`, `mon.names`, and `pos.*` indices
    automatically, eliminating a common source of bugs.

4.  **Prefer explicit parameterization.** Use `NormalV` for variance
    parameterization and `NormalP` for precision, rather than relying on
    implicit conventions. The DSL enforces this explicitness.

5.  **Leverage `initial_values()` as a starting point.** The generated
    initial values are reasonable defaults. For difficult posterior
    geometries, you may need to adjust them manually.

6.  **Check the DAG.** `plot(spec)` provides a quick visual sanity check
    that dependencies are correctly specified.

7.  **For complex models, build incrementally.** Start with a simple
    version of the model that compiles, verify it runs, then add
    complexity one statement at a time.

## Supported algorithms

The compiled Model function is compatible with all lucifer estimation
algorithms: iterative quadrature
([`IterativeQuadrature()`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md)),
Laplace approximation
([`LaplaceApproximation()`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md)),
all 55 MCMC algorithms (via
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md)),
Population Monte Carlo
([`PMC()`](https://robustecologies.github.io/lucifer/reference/PMC.md)),
and Variational Bayes
([`VariationalBayes()`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md)).
No modification is needed; the generated function follows the universal
Model function interface.

``` r

# NUTS (default recommendation)
fit_nuts <- lucifer(spec$Model, Data, IV, Iterations = 10000, Algorithm = "NUTS")

# Adaptive Metropolis
fit_am <- lucifer(spec$Model, Data, IV, Iterations = 50000, Algorithm = "AM")

# Laplace approximation
fit_la <- LaplaceApproximation(spec$Model, Data, parm = IV, Method = "SPG")

# Variational Bayes
fit_vb <- VariationalBayes(spec$Model, Data, parm = IV, Method = "Salimans1")
```

## Technical reference

### Generated code structure

The compiler emits code blocks in this fixed order, matching lucifer
convention:

1.  **Parameter extraction**: `param <- parm[Data$pos.param]`
2.  **Constraint application**:
    `param <- interval(parm[Data$pos.param], lower, upper); parm[Data$pos.param] <- param`
3.  **Log-priors** (hyperpriors first, then priors):
    `param.prior <- sum(dfun(param, args, log = TRUE))`
4.  **Deterministic relationships**:
    `lhs <- expression_with_Data_prefix`
5.  **Log-likelihood**:
    `LL <- sum(dfun(Data$response, args, log = TRUE))`
6.  **Log-posterior assembly**: `LP <- LL + prior1 + prior2 + ...`
7.  **Posterior predictive**: `yhat <- rfun(N, args)`
8.  **Return list**:
    `list(LP = LP, Dev = -2 * LL, Monitor = LP, yhat = yhat, parm = parm)`

### Data\$pos.\* convention

Position indices are stored in the Data list as `Data$pos.param_name`.
For scalar parameters, this is a single integer. For vectorized
parameters (e.g., `beta[j]` with `j = 1,...,J`), this is an integer
vector `1:J`. These indices map directly into the `parm` vector.

### Intermediate representation

The IR (class `model_ir`) is a list of node objects, each containing:

``` r

# Stochastic node
list(
    type = "stochastic",
    lhs = "beta",
    lhs_subscript = "j",
    index = list(var = "j", range = "1:J"),
    distribution = "NormalV",
    args = c("0", "1000"),
    role = "prior",          # "likelihood" | "prior" | "hyperprior"
    constraint = list(type = "none", lower = -Inf, upper = Inf)
)

# Deterministic node
list(
    type = "deterministic",
    lhs = "mu",
    expression = "X %*% beta",
    refs = c("X", "beta")
)
```

The IR is accessible via `spec$ir` and is useful for programmatic model
manipulation or debugging.

## Session info

``` r

sessionInfo()
```
