# Working with compiled C++ models

## Introduction

Every MCMC iteration in lucifer requires evaluating the log-posterior
and, for gradient-based samplers like NUTS, its gradient. When the model
is an R function, each evaluation crosses the R/C++ boundary: the C++
sampler marshals parameters into R objects, calls the R function, and
unmarshals the result. For NUTS with \\n\\ parameters, finite-difference
gradients require \\2n + 1\\ such round-trips per gradient evaluation,
and a single iteration with \\L\\ leapfrog steps demands \\L \times
(2n + 1)\\ R function calls. With 10 parameters and \\L = 30\\, that is
over 600 R calls per iteration.

The compiled model backend eliminates this bottleneck. When a model is
specified via
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md),
the system can generate a self-contained C++ source file that evaluates
the log-posterior and its analytical gradient without any R callbacks.
The generated code is compiled at runtime via
[`Rcpp::sourceCpp()`](https://rdrr.io/pkg/Rcpp/man/sourceCpp.html),
cached for reuse, and integrated transparently with lucifer’s sampler
infrastructure. The result: gradient computation becomes a single C++
function call instead of hundreds of R round-trips, yielding speedups of
100–500x for gradient evaluation and 5–15x for end-to-end NUTS sampling.

This vignette demonstrates the complete workflow, validates correctness
against R models, benchmarks performance across model types and sizes,
and documents the full API.

## Quick start

The simplest path from model specification to compiled C++ sampling
requires three steps: specify the model with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md),
compile it with
[`compile()`](https://robustecologies.github.io/lucifer/reference/compile.md),
and pass the compiled object to
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

``` r

# 1. Specify model
spec <- model_spec("
    y ~ Normal(mu, sigma)
    mu = X %*% beta
    beta[j] ~ Normal(0, 10), j = 1,...,J
    sigma ~ HalfCauchy(2.5)
", backend = "cpp")

# 2. Generate data
set.seed(42)
N <- 100; J <- 3
X <- cbind(1, matrix(rnorm(N * (J - 1)), N, J - 1))
true_beta <- c(2, -1, 0.5)
true_sigma <- 1.5
y <- as.numeric(X %*% true_beta + rnorm(N, 0, true_sigma))

# 3. Compile to C++
Data <- spec$data_template(y = y, X = X, J = J)
compiled <- compile(spec, Data, verbose = FALSE)
print(compiled)
```

The compiled object can be passed directly to
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
which automatically extracts the analytical gradient and registers
native function pointers for the sampler loop.

``` r

iv <- compiled$initial_values(Data)
fit <- lucifer(compiled, Data, iv,
               Iterations = 2000, Status = 2000,
               Thinning = 1, Chains = 1,
               Algorithm = "NUTS",
               Specs = list(A = 1000, delta = 0.65,
                            epsilon = 0.1, Lmax = 10))
```

Since the compiled model operates in unconstrained space
(log-transformed for positive parameters), posterior samples need
back-transformation for interpretation.

``` r

post <- fit$Posterior2
post[, J + 1] <- exp(post[, J + 1])  # sigma = exp(log_sigma)
colnames(post) <- c(paste0("beta[", 1:J, "]"), "sigma")

knitr::kable(
    data.frame(
        parameter = colnames(post),
        true = c(true_beta, true_sigma),
        mean = round(colMeans(post), 3),
        sd = round(apply(post, 2, sd), 3),
        q025 = round(apply(post, 2, quantile, 0.025), 3),
        q975 = round(apply(post, 2, quantile, 0.975), 3)
    ),
    caption = "Posterior summary from compiled NUTS (2000 iterations)",
    row.names = FALSE
)
```

``` r

par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
for (j in 1:ncol(post)) {
    plot(post[, j], type = "l", col = "#2E5090",
         ylab = colnames(post)[j], xlab = "Iteration",
         main = colnames(post)[j])
    abline(h = c(true_beta, true_sigma)[j], col = "#E07020", lwd = 2, lty = 2)
}
```

## How it works

The compilation pipeline transforms the model specification through four
stages, each producing a well-defined intermediate product.

**Stage 1: parsing and IR construction.** The
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
function parses the probabilistic notation into a topologically sorted
directed acyclic graph (DAG). Each node is classified as data,
parameter, hyperparameter, or deterministic quantity. Constraint types
(positive, unit, bounded) are inferred from distribution support. This
stage is identical for R and C++ backends.

**Stage 2: C++ code generation.** The
[`compile()`](https://robustecologies.github.io/lucifer/reference/compile.md)
function walks the DAG and emits a self-contained C++ source file
containing three functions: `init_data_HASH()` populates a static data
struct from the R Data list; `eval_model_HASH()` computes the
log-posterior with Jacobian-adjusted parameter transforms;
`grad_model_HASH()` computes the analytical gradient via
distribution-specific formulas and chain rule through deterministic
compositions.

**Stage 3: runtime compilation.** The generated source is compiled via
[`Rcpp::sourceCpp()`](https://rdrr.io/pkg/Rcpp/man/sourceCpp.html) and
cached in `tools::R_user_dir("lucifer", "cache")`. Subsequent calls with
the same model specification skip compilation entirely.

**Stage 4: native dispatch.** When
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
receives a compiled model, it registers C-style function pointers in the
sampler infrastructure. The `call_model()` and `compute_gradient()`
functions in the C++ sampler loop detect these pointers and dispatch
directly, bypassing the R function call mechanism, marshaling, and RNG
synchronization overhead.

### Inspecting the generated code

The
[`code()`](https://robustecologies.github.io/lucifer/reference/code.md)
generic method provides structured access to the generated C++ source.
The `section` argument selects which part to display:

``` r

code(compiled, section = "eval")
```

The gradient function mirrors the evaluation structure, with the chain
rule applied in reverse through the DAG:

``` r

code(compiled, section = "grad")
```

The full source can be saved to disk for inspection or manual editing
via `code(compiled, file = "my_model.cpp")`.

## Unconstrained parameterization

A critical design decision: the compiled model operates entirely in
unconstrained space. Parameters with positivity constraints are
log-transformed; unit-interval parameters use the logit transform;
bounded parameters use the scaled logit. The log-Jacobian determinant of
the transformation is added to the log-posterior so that HMC/NUTS
samples the correct target density.

This differs from the R
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
backend, which uses
[`interval()`](https://robustecologies.github.io/lucifer/reference/interval.md)
(clamping) to enforce constraints. Clamping has zero gradient at
boundaries and introduces discontinuities that break gradient-based
samplers. The differentiable bijections used by the C++ backend are
essential for correct HMC dynamics.

The transforms and their Jacobians, implemented in
`inst/include/lucifer/transforms.h`, are:

\\ \begin{aligned} \text{Positive: } & x = e^y, \quad \log\|J\| = y \\
\text{Unit: } & x = \text{logit}^{-1}(y), \quad \log\|J\| = \log x +
\log(1-x) \\ \text{Bounded } \[a,b\]\text{: } & x = a +
(b-a)\\\text{logit}^{-1}(y), \quad \log\|J\| = \log(b-a) + \log x' +
\log(1-x') \end{aligned} \\

where \\x' = (x-a)/(b-a)\\. The `initial_values()` function of a
compiled model automatically converts constrained starting values to
unconstrained space using the inverse transforms.

## Model examples

### Normal linear regression

The quickstart example above demonstrates linear regression. Here we
verify gradient correctness by comparing the analytical gradient against
central finite differences.

``` r

iv <- compiled$initial_values(Data)
analytic <- compiled$Grad(iv, Data)
eps <- 1e-5
fd <- numeric(length(iv))
for (i in seq_along(iv)) {
    p <- m <- iv
    p[i] <- iv[i] + eps; m[i] <- iv[i] - eps
    fd[i] <- (compiled$Model(p, Data)$LP - compiled$Model(m, Data)$LP) / (2 * eps)
}

knitr::kable(
    data.frame(
        parameter = c(paste0("beta_unc[", 1:J, "]"), "log_sigma"),
        analytic = analytic,
        finite_diff = fd,
        abs_error = abs(analytic - fd)
    ),
    digits = 8, row.names = FALSE,
    caption = "Gradient validation: analytical vs central finite differences"
)
```

### Poisson GLM with log link

The Poisson GLM uses `lambda = exp(X %*% beta)`, which requires the
gradient to chain-rule through the exponential function. The compiled
backend handles this automatically: the derivative of `exp(z)` evaluated
at `z = X %*% beta` is `exp(z) = \lambda` itself, so the gradient
contribution from the likelihood is \\X^\top \cdot (y - \lambda)\\ per
observation.

``` r

spec_pois <- model_spec("
    y ~ Poisson(lambda)
    lambda = exp(X %*% beta)
    beta[j] ~ Normal(0, 5), j = 1,...,J
")

set.seed(10)
N_p <- 80; J_p <- 3
X_p <- cbind(1, matrix(rnorm(N_p * (J_p - 1)), N_p, J_p - 1))
true_beta_p <- c(0.5, 0.3, -0.2)
y_p <- rpois(N_p, exp(X_p %*% true_beta_p))

Data_p <- spec_pois$data_template(y = y_p, X = X_p, J = J_p)
compiled_p <- compile(spec_pois, Data_p, verbose = FALSE)
iv_p <- compiled_p$initial_values(Data_p)

# Gradient check
grad_p <- compiled_p$Grad(iv_p, Data_p)
fd_p <- numeric(length(iv_p))
for (i in seq_along(iv_p)) {
    p <- m <- iv_p; p[i] <- iv_p[i] + eps; m[i] <- iv_p[i] - eps
    fd_p[i] <- (compiled_p$Model(p, Data_p)$LP -
                compiled_p$Model(m, Data_p)$LP) / (2 * eps)
}

knitr::kable(
    data.frame(
        parameter = paste0("beta[", seq_along(iv_p), "]"),
        analytic = grad_p,
        finite_diff = fd_p,
        abs_error = abs(grad_p - fd_p)
    ),
    digits = 6, row.names = FALSE,
    caption = "Poisson GLM gradient validation (chain rule through exp)"
)
```

``` r

fit_p <- lucifer(compiled_p, Data_p, iv_p,
                 Iterations = 2000, Status = 2000,
                 Thinning = 1, Chains = 1,
                 Algorithm = "NUTS",
                 Specs = list(A = 1000, delta = 0.65,
                              epsilon = 0.05, Lmax = 10))
post_p <- fit_p$Posterior2

knitr::kable(
    data.frame(
        parameter = paste0("beta[", 1:J_p, "]"),
        true = true_beta_p,
        posterior_mean = round(colMeans(post_p), 3),
        posterior_sd = round(apply(post_p, 2, sd), 3)
    ),
    row.names = FALSE,
    caption = paste0("Poisson GLM posterior (acceptance: ",
                     round(fit_p$Acceptance.Rate, 3), ")")
)
```

### Logistic regression with invlogit link

Binary response models use `p = invlogit(X %*% beta)`, requiring the
chain rule through the inverse logit function. The derivative is
\\\text{invlogit}'(z) = p(1-p)\\, which the gradient generator applies
automatically.

``` r

spec_logit <- model_spec("
    y ~ Bernoulli(p)
    p = invlogit(X %*% beta)
    beta[j] ~ Normal(0, 5), j = 1,...,J
")

set.seed(20)
N_l <- 120; J_l <- 3
X_l <- cbind(1, matrix(rnorm(N_l * (J_l - 1)), N_l, J_l - 1))
true_beta_l <- c(-0.5, 1.0, -0.3)
p_true <- plogis(X_l %*% true_beta_l)
y_l <- rbinom(N_l, 1, p_true)

Data_l <- spec_logit$data_template(y = y_l, X = X_l, J = J_l)
compiled_l <- compile(spec_logit, Data_l, verbose = FALSE)
iv_l <- compiled_l$initial_values(Data_l)

# Gradient check
grad_l <- compiled_l$Grad(iv_l, Data_l)
fd_l <- numeric(length(iv_l))
for (i in seq_along(iv_l)) {
    p <- m <- iv_l; p[i] <- iv_l[i] + eps; m[i] <- iv_l[i] - eps
    fd_l[i] <- (compiled_l$Model(p, Data_l)$LP -
                compiled_l$Model(m, Data_l)$LP) / (2 * eps)
}

knitr::kable(
    data.frame(
        parameter = paste0("beta[", seq_along(iv_l), "]"),
        analytic = grad_l, finite_diff = fd_l,
        abs_error = abs(grad_l - fd_l)
    ),
    digits = 6, row.names = FALSE,
    caption = "Logistic regression gradient validation (chain rule through invlogit)"
)
```

``` r

fit_l <- lucifer(compiled_l, Data_l, iv_l,
                 Iterations = 2000, Status = 2000,
                 Thinning = 1, Chains = 1,
                 Algorithm = "NUTS",
                 Specs = list(A = 1000, delta = 0.65,
                              epsilon = 0.05, Lmax = 10))
post_l <- fit_l$Posterior2

knitr::kable(
    data.frame(
        parameter = paste0("beta[", 1:J_l, "]"),
        true = true_beta_l,
        posterior_mean = round(colMeans(post_l), 3),
        posterior_sd = round(apply(post_l, 2, sd), 3)
    ),
    row.names = FALSE,
    caption = paste0("Logistic regression posterior (acceptance: ",
                     round(fit_l$Acceptance.Rate, 3), ")")
)
```

### Gamma regression

The Gamma likelihood demonstrates the system’s ability to handle
distributions with two positive parameters, both requiring log
transforms.

``` r

spec_gam <- model_spec("
    y ~ Gamma(shape, rate)
    shape ~ Exponential(1)
    rate ~ Exponential(1)
")

set.seed(30)
true_shape <- 2.5; true_rate <- 1.2
y_g <- rgamma(60, shape = true_shape, rate = true_rate)
Data_g <- spec_gam$data_template(y = y_g)
compiled_g <- compile(spec_gam, Data_g, verbose = FALSE)
iv_g <- compiled_g$initial_values(Data_g)

fit_g <- lucifer(compiled_g, Data_g, iv_g,
                 Iterations = 2000, Status = 2000,
                 Thinning = 1, Chains = 1,
                 Algorithm = "NUTS",
                 Specs = list(A = 1000, delta = 0.65,
                              epsilon = 0.05, Lmax = 10))
post_g <- exp(fit_g$Posterior2)  # Both params positive

knitr::kable(
    data.frame(
        parameter = c("shape", "rate"),
        true = c(true_shape, true_rate),
        posterior_mean = round(colMeans(post_g), 3),
        posterior_sd = round(apply(post_g, 2, sd), 3)
    ),
    row.names = FALSE,
    caption = paste0("Gamma model posterior (acceptance: ",
                     round(fit_g$Acceptance.Rate, 3), ")")
)
```

### Beta regression

Unit-interval parameters use the logit transform, which is the most
numerically delicate bijection. The compiled backend uses a numerically
stable sigmoid implementation that handles extreme values of the
unconstrained parameter without overflow.

``` r

spec_beta <- model_spec("
    y ~ Beta(a, b)
    a ~ Gamma(2, 1)
    b ~ Gamma(2, 1)
")

set.seed(40)
true_a <- 3; true_b <- 7
y_b <- rbeta(50, shape1 = true_a, shape2 = true_b)
Data_b <- spec_beta$data_template(y = y_b)
compiled_b <- compile(spec_beta, Data_b, verbose = FALSE)
iv_b <- compiled_b$initial_values(Data_b)

fit_b <- lucifer(compiled_b, Data_b, iv_b,
                 Iterations = 2000, Status = 2000,
                 Thinning = 1, Chains = 1,
                 Algorithm = "NUTS",
                 Specs = list(A = 1000, delta = 0.65,
                              epsilon = 0.05, Lmax = 10))
post_b <- exp(fit_b$Posterior2)

knitr::kable(
    data.frame(
        parameter = c("a", "b"),
        true = c(true_a, true_b),
        posterior_mean = round(colMeans(post_b), 3),
        posterior_sd = round(apply(post_b, 2, sd), 3)
    ),
    row.names = FALSE,
    caption = paste0("Beta model posterior (acceptance: ",
                     round(fit_b$Acceptance.Rate, 3), ")")
)
```

## Comprehensive benchmarks

The performance advantage of the compiled backend comes from two
sources: compiled model evaluation (eliminating R marshaling) and
analytical gradients (eliminating finite-difference round-trips). The
following benchmarks measure each component separately and their
combined effect on NUTS sampling.

### Model evaluation benchmark

We compare the wall-clock time for repeated model evaluations across
increasing data sizes, holding the number of parameters fixed at \\J + 1
= 6\\. Each configuration is measured as the median of 3 timing runs,
each with 5000 evaluations, to reduce noise from OS scheduling and
garbage collection.

``` r

spec_bench <- model_spec("
    y ~ Normal(mu, sigma)
    mu = X %*% beta
    beta[j] ~ Normal(0, 10), j = 1,...,J
    sigma ~ HalfCauchy(2.5)
")

J_bench <- 5
N_values <- c(50, 100, 200, 500, 1000)

# Helper: median of k timing runs for stability.
# fn must be a zero-argument function (closure).
bench_median <- function(fn, n_reps, k = 3) {
    times <- numeric(k)
    for (j in seq_len(k))
        times[j] <- system.time(for (i in seq_len(n_reps)) fn())[3]
    median(times) / n_reps
}

set.seed(666)
eval_results <- do.call(rbind, lapply(N_values, function(N) {
    X <- cbind(1, matrix(rnorm(N * (J_bench - 1)), N, J_bench - 1))
    yy <- as.numeric(X %*% rnorm(J_bench) + rnorm(N))
    D <- spec_bench$data_template(y = yy, X = X, J = J_bench)
    cc <- compile(spec_bench, D, verbose = FALSE)
    iv_r <- spec_bench$initial_values(D)
    iv_c <- cc$initial_values(D)

    # Warmup
    for (w in 1:50) { spec_bench$Model(iv_r, D); cc$Model(iv_c, D) }

    t_r <- bench_median(function() spec_bench$Model(iv_r, D), n_reps = 5000)
    t_c <- bench_median(function() cc$Model(iv_c, D), n_reps = 5000)
    data.frame(N = N, R_ms = t_r * 1000,
               Cpp_ms = t_c * 1000,
               speedup = t_r / max(t_c, 1e-9))
}))

knitr::kable(eval_results, digits = c(0, 4, 4, 1), row.names = FALSE,
             col.names = c("N", "R (ms/call)", "C++ (ms/call)", "Speedup"),
             caption = "Model evaluation: median of 3 runs, 5000 calls each (J = 5)")
```

``` r

eval_long <- rbind(
    data.frame(N = eval_results$N, backend = "R model", ms = eval_results$R_ms),
    data.frame(N = eval_results$N, backend = "C++ compiled", ms = eval_results$Cpp_ms)
)
ggplot(eval_long, aes(x = N, y = ms, color = backend)) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    scale_color_manual(values = c("R model" = "#2E5090",
                                   "C++ compiled" = "#E07020")) +
    labs(x = "Number of observations (N)", y = "Time per evaluation (ms)",
         color = "Backend",
         title = "Model evaluation cost",
         subtitle = paste0("J = ", J_bench, " parameters, median of 3 runs")) +
    theme_relab() +
    theme(legend.position = "bottom")
```

The C++ compiled model is consistently faster, with the advantage
growing as \\N\\ increases because the compiled code eliminates both the
R boundary crossing overhead and the interpretive overhead of R-level
distribution calls.

### Gradient computation benchmark

The gradient speedup is the dominant performance advantage. The R model
requires \\2n + 1\\ model evaluations for central finite differences;
the compiled model computes the exact analytical gradient in a single
C++ call. We fix \\N = 200\\ and increase the number of parameters \\J\\
to show the linear scaling. The theoretical speedup is \\2n + 1\\ (where
\\n\\ is the number of parameters), since finite differences require
exactly that many model evaluations while the analytical gradient
requires one.

``` r

N_grad <- 200
J_values <- c(3, 5, 8, 12, 20)

set.seed(666)
grad_results <- do.call(rbind, lapply(J_values, function(J) {
    X <- cbind(1, matrix(rnorm(N_grad * (J - 1)), N_grad, J - 1))
    yy <- as.numeric(X %*% rnorm(J) + rnorm(N_grad))
    D <- spec_bench$data_template(y = yy, X = X, J = J)
    cc <- compile(spec_bench, D, verbose = FALSE)
    iv_r <- spec_bench$initial_values(D)
    iv_c <- cc$initial_values(D)
    n_parm <- J + 1

    # Warmup
    for (w in 1:20) { spec_bench$Model(iv_r, D); cc$Grad(iv_c, D) }

    # Finite-difference gradient (R)
    n_fd <- max(100, 500 %/% n_parm)
    fd_fn <- function() {
        g <- numeric(length(iv_r))
        for (k in seq_along(iv_r)) {
            p <- m <- iv_r
            p[k] <- iv_r[k] + 1e-6; m[k] <- iv_r[k] - 1e-6
            g[k] <- (spec_bench$Model(p, D)$LP -
                     spec_bench$Model(m, D)$LP) / 2e-6
        }
    }
    t_fd <- bench_median(fd_fn, n_reps = n_fd) * 1000

    # Analytical gradient (C++): many reps since each call is ~5us
    t_cpp <- bench_median(function() cc$Grad(iv_c, D), n_reps = 20000) * 1000

    data.frame(J = J, n_parm = n_parm,
               FD_ms = t_fd, Cpp_ms = t_cpp,
               speedup = t_fd / max(t_cpp, 1e-6))
}))

knitr::kable(grad_results, digits = c(0, 0, 3, 4, 0), row.names = FALSE,
             col.names = c("J", "n_parm", "FD (ms/call)", "C++ (ms/call)", "Speedup"),
             caption = "Gradient computation per call (N = 200, median of 3 runs)")
```

``` r

grad_long <- rbind(
    data.frame(n_parm = grad_results$n_parm,
               method = "R finite differences", ms = grad_results$FD_ms),
    data.frame(n_parm = grad_results$n_parm,
               method = "C++ analytical", ms = grad_results$Cpp_ms)
)
ggplot(grad_long, aes(x = n_parm, y = ms, color = method)) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    scale_y_log10() +
    scale_color_manual(values = c("R finite differences" = "#2E5090",
                                   "C++ analytical" = "#E07020")) +
    labs(x = "Number of parameters",
         y = "Time per gradient call (ms, log scale)",
         color = "Method",
         title = "Gradient computation cost",
         subtitle = "N = 200 observations, median of 3 timing runs") +
    theme_relab() +
    theme(legend.position = "bottom")
```

Finite-difference cost grows linearly with the number of parameters
(each parameter requires two model evaluations), while the analytical
gradient cost remains nearly constant. The theoretical speedup therefore
scales approximately as \\2n + 1\\. The ratio to theoretical is in
practice ~7x because the theoretical formula only counts model
evaluations, while in practice each R model evaluation has more overhead
than the C++ gradient call (R interpreter, marshaling, PutRNGstate,
etc.).:

``` r

# Add theoretical speedup line
grad_results$theoretical <- 2 * grad_results$n_parm + 1

speedup_long <- rbind(
    data.frame(n_parm = grad_results$n_parm,
               type = "Measured", speedup = grad_results$speedup),
    data.frame(n_parm = grad_results$n_parm,
               type = "Theoretical (2n+1)", speedup = grad_results$theoretical)
)

ggplot(speedup_long, aes(x = n_parm, y = speedup, color = type,
                          linetype = type)) +
    geom_line(linewidth = 1.1) +
    geom_point(data = speedup_long[speedup_long$type == "Measured", ],
               size = 3.5) +
    geom_text(data = speedup_long[speedup_long$type == "Measured", ],
              aes(label = paste0(round(speedup), "x")),
              vjust = -1.2, size = 3.5, fontface = "bold", show.legend = FALSE) +
    scale_color_manual(values = c("Measured" = "#C04040",
                                   "Theoretical (2n+1)" = "grey50")) +
    scale_linetype_manual(values = c("Measured" = "solid",
                                      "Theoretical (2n+1)" = "dashed")) +
    labs(x = "Number of parameters",
         y = "Speedup (FD time / analytical time)",
         color = NULL, linetype = NULL,
         title = "Gradient speedup: analytical C++ vs finite differences",
         subtitle = "Dashed line: theoretical 2n+1 scaling") +
    theme_relab() +
    expand_limits(y = max(grad_results$speedup) * 1.1)
```

### End-to-end NUTS benchmark

The benchmark that matters most is end-to-end sampling time with NUTS.
This combines model evaluation, gradient computation, leapfrog
integration, and tree-building.

``` r

nuts_results <- do.call(rbind, lapply(N_values[1:4], function(N) {
    set.seed(666)
    X <- cbind(1, matrix(rnorm(N * (J_bench - 1)), N, J_bench - 1))
    yy <- as.numeric(X %*% rnorm(J_bench) + rnorm(N))
    D <- spec_bench$data_template(y = yy, X = X, J = J_bench)
    cc <- compile(spec_bench, D, verbose = FALSE)
    iv_r <- spec_bench$initial_values(D)
    iv_c <- cc$initial_values(D)
    iters <- 2000
    nuts_specs <- list(A = 1000, delta = 0.65, epsilon = 0.1, Lmax = 10)

    t_r <- system.time({
        fit_r <- lucifer(spec_bench, D, iv_r,
                         Iterations = iters, Status = iters,
                         Thinning = 1, Chains = 1,
                         Algorithm = "NUTS", Specs = nuts_specs)
    })[3]

    t_c <- system.time({
        fit_c <- lucifer(cc, D, iv_c,
                         Iterations = iters, Status = iters,
                         Thinning = 1, Chains = 1,
                         Algorithm = "NUTS", Specs = nuts_specs)
    })[3]

    data.frame(N = N, R_s = t_r, Cpp_s = t_c,
               speedup = t_r / max(t_c, 0.01),
               R_accept = fit_r$Acceptance.Rate,
               Cpp_accept = fit_c$Acceptance.Rate)
}))

knitr::kable(nuts_results, digits = c(0, 2, 2, 1, 3, 3), row.names = FALSE,
             col.names = c("N", "R (s)", "C++ (s)", "Speedup",
                           "R accept", "C++ accept"),
             caption = "NUTS end-to-end wall-clock time (2000 iterations, J = 5)")
```

``` r

nuts_long <- rbind(
    data.frame(N = nuts_results$N, backend = "R (finite diff)",
               seconds = nuts_results$R_s),
    data.frame(N = nuts_results$N, backend = "C++ (analytical)",
               seconds = nuts_results$Cpp_s)
)

ggplot(nuts_long, aes(x = factor(N), y = seconds, fill = backend)) +
    geom_col(position = "dodge", width = 0.6) +
    scale_fill_manual(values = c("R (finite diff)" = "#2E5090",
                                  "C++ (analytical)" = "#E07020")) +
    labs(x = "Number of observations (N)", y = "Wall-clock time (seconds)",
         fill = "Backend",
         title = "NUTS sampling time: R vs compiled C++",
         subtitle = "2000 iterations, J = 5 parameters, single chain") +
    theme_relab() +
    theme(legend.position = "bottom")
```

## Parallel chains

Compiled models support multi-chain parallel sampling. When
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
receives a compiled model with `Chains > 1`, it extracts the
serializable metadata (C++ source code, model hash, IR), passes it to
each child process, and recompiles from the cache in each child. Since
[`Rcpp::sourceCpp()`](https://rdrr.io/pkg/Rcpp/man/sourceCpp.html)
detects the cached compiled object, recompilation takes approximately
100ms per child rather than the full 3–5 seconds, making the overhead
negligible relative to sampling time.

``` r

set.seed(42)
fit_par <- lucifer(compiled, Data, iv,
                   Iterations = 1500, 
                   Thinning = 1, Chains = 3, CPUs = 3,
                   Algorithm = "NUTS",
                   Specs = list(A = 750, delta = 0.65,
                                epsilon = 0.1, Lmax = 10))
```

``` r

knitr::kable(
    data.frame(
        chain = 1:3,
        acceptance = round(sapply(fit_par$Chains,
            function(c) c$Acceptance.Rate), 3)
    ),
    row.names = FALSE,
    caption = "Per-chain acceptance rates (compiled NUTS, 3 chains)"
)
```

The combined posterior from all three chains is available in
`fit_par$Posterior2`. Back-transforming the constrained parameters and
comparing against truth:

``` r

post_par <- fit_par$Posterior2
post_par[, J + 1] <- exp(post_par[, J + 1])
colnames(post_par) <- c(paste0("beta[", 1:J, "]"), "sigma")

knitr::kable(
    data.frame(
        parameter = colnames(post_par),
        true = c(true_beta, true_sigma),
        mean = round(colMeans(post_par), 3),
        sd = round(apply(post_par, 2, sd), 3)
    ),
    row.names = FALSE,
    caption = "Combined posterior from 3 parallel compiled chains"
)
```

## NUTS diagnostics

The NUTS sampler tracks per-iteration diagnostics that are critical for
assessing sampling quality: Hamiltonian energy, tree depth, divergent
transitions, and leapfrog step counts. These are stored in the
`$NUTS.Diagnostics` field of the fitted object and can be inspected with
[`check_nuts()`](https://robustecologies.github.io/lucifer/reference/check_nuts.md),
which provides a summary similar to Stan’s
[`check_hmc_diagnostics()`](https://mc-stan.org/rstan/reference/check_hmc_diagnostics.html).

``` r

spec_diag <- model_spec("
    y ~ Normal(mu, sigma)
    mu = X %*% beta
    beta[j] ~ Normal(0, 10), j = 1,...,J
    sigma ~ HalfCauchy(2.5)
")

set.seed(666)
N_d <- 100; J_d <- 4
X_d <- cbind(1, matrix(rnorm(N_d * (J_d - 1)), N_d, J_d - 1))
y_d <- as.numeric(X_d %*% c(1, -0.5, 0.3, 0.8) + rnorm(N_d, 0, 2))
Data_d <- spec_diag$data_template(y = y_d, X = X_d, J = J_d)
compiled_d <- compile(spec_diag, Data_d, verbose = FALSE)
iv_d <- compiled_d$initial_values(Data_d)

fit_d <- lucifer(compiled_d, Data_d, iv_d,
                 Iterations = 2000, Status = 2000,
                 Thinning = 1, Chains = 3, CPUs = 3,
                 Algorithm = "NUTS",
                 Specs = list(A = 1000, delta = 0.8,
                              epsilon = 0.1, Lmax = 12))

nuts_diag <- check_nuts(fit_d, warmup = 200)
nuts_diag
```

The
[`check_nuts()`](https://robustecologies.github.io/lucifer/reference/check_nuts.md)
function returns an S3 `nuts_check` object. When the fit contains
multiple chains, diagnostics are reported per chain (E-BFMI) or
aggregated with chain structure preserved for visualization. The
`warmup` argument trims initial adaptation iterations from the
diagnostic plots (but not from the summary statistics), removing the
transient spikes that occur during step-size adaptation. The
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) method produces
a 4-panel diagnostic summary where each chain is rendered in a distinct
color, making it easy to assess convergence and mixing visually.
Divergent transitions are highlighted in red across all panels.

``` r

plot(nuts_diag)
```

## Auto-compilation via `backend`

When
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
is called with `backend = "cpp"` or `backend = "auto"`, the
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
function automatically compiles the model the first time it is called
with Data. This makes the compiled backend fully transparent: the user
writes the same notation and the system selects the fastest available
implementation.

``` r

spec_auto <- model_spec("
    y ~ Normal(mu, sigma)
    mu ~ Normal(0, 100)
    sigma ~ HalfCauchy(25)
", backend = "auto")

set.seed(99)
Data_auto <- spec_auto$data_template(y = rnorm(30, 5, 2))
iv_auto <- spec_auto$initial_values(Data_auto)

fit_auto <- lucifer(spec_auto, Data_auto, iv_auto,
                    Iterations = 1000, Status = 1000,
                    Thinning = 1, Chains = 1,
                    Algorithm = "NUTS",
                    Specs = list(A = 500, delta = 0.65,
                                 epsilon = 0.1, Lmax = 10))
```

The `backend = "auto"` setting attempts C++ compilation silently and
falls back to the R backend if compilation fails (missing compiler,
unsupported distribution, etc.). No user intervention is needed.

## Supported distributions

The C++ gradient library provides analytical gradients for 24
distributions, covering the vast majority of model_spec usage. Each
distribution has gradients w.r.t. all arguments (the random variable and
all distribution parameters), enabling chain-rule propagation through
the model DAG.

``` r

supported <- data.frame(
    distribution = c("Normal", "NormalV", "NormalP", "StudentT",
                     "LogNormal", "LogNormalP", "Laplace", "LaplaceP",
                     "Weibull",
                     "HalfCauchy", "HalfNormal", "HalfT",
                     "Gamma", "InvGamma", "Exponential", "Pareto",
                     "InvChiSq", "InvGaussian",
                     "Beta", "Uniform",
                     "Poisson", "Bernoulli", "Binomial", "NegBinomial"),
    constraint = c(rep("none", 9),
                   rep("positive", 9),
                   "unit", "bounded",
                   rep("none (discrete)", 4)),
    category = c(rep("Continuous", 9),
                 rep("Positive continuous", 9),
                 "Bounded continuous", "Bounded continuous",
                 rep("Discrete", 4))
)
knitr::kable(supported, row.names = FALSE,
             caption = "Distributions with C++ analytical gradients (24 of 86 total)")
```

Models using distributions not in this list (Stable, MVN, Wishart,
copulas, etc.) will fall back to the R backend with a clear diagnostic
message. The gradient header also includes `MVNGrad` and `DirichletGrad`
structs for multivariate distributions; wiring these into the code
generator for compiled multivariate models is planned for a future
release.

## Chain rule through compositions

The gradient generator handles deterministic nodes that wrap matrix
multiplication in common link functions. The supported compositions are:

- **Identity**: `mu = X %*% beta` \\\to\\ \\\partial/\partial\beta =
  X^\top \cdot (\partial\text{LL}/\partial\mu)\\
- **Exponential**: `lambda = exp(X %*% beta)` \\\to\\ multiplies by
  \\\lambda\\ (Poisson, log-link GLMs)
- **Inverse logit**: `p = invlogit(X %*% beta)` \\\to\\ multiplies by
  \\p(1-p)\\ (logistic regression)
- **Logarithm**: `z = log(X %*% beta)` \\\to\\ divides by argument
- **Square root**: `s = sqrt(X %*% beta)` \\\to\\ multiplies by
  \\1/(2s)\\

Deeper compositions (e.g., `log(1 + exp(x))`, `sin(2*pi*t/12)`) are not
yet handled by the symbolic gradient system but are supported by the
dual number automatic differentiation library (see below). Full
integration of the AD library into the code generator is planned for a
future release.

## Automatic differentiation library

The package includes a header-only forward-mode automatic
differentiation library at `inst/include/lucifer/dual.h`. It implements
dual numbers \\(v, \dot{v})\\ where \\v\\ is the function value and
\\\dot{v}\\ is the derivative w.r.t. a seed variable. All arithmetic
operators and 20+ elementary functions propagate derivatives
automatically via the chain rule.

The library supports: `exp`, `log`, `log1p`, `sqrt`, `pow`, `abs`,
`sin`, `cos`, `tan`, `tanh`, `sinh`, `cosh`, `asin`, `acos`, `atan`,
`invlogit`, `logit`, and `softplus`. It is designed for integration into
generated model code when the symbolic gradient system cannot handle a
deterministic expression.

As an example, the derivative of the softplus function \\\log(1 +
e^{ax+b})\\ is computed automatically:

``` r

# The softplus derivative is a*sigmoid(a*x + b)
x_ad <- 1.5; a_ad <- 2.0; b_ad <- -1.0
r_ad <- lucifer:::test_dual_softplus(x_ad, a_ad, b_ad)
expected_ad <- a_ad * plogis(a_ad * x_ad + b_ad)

# Seasonal derivative: d/dx sin(2*pi*x/12) = 2*pi/12 * cos(2*pi*x/12)
x_seas <- 3.0
r_seas <- lucifer:::test_dual_seasonal(x_seas)
expected_seas <- 2 * pi / 12 * cos(2 * pi * x_seas / 12)

knitr::kable(
    data.frame(
        expression = c("log(1 + exp(2x - 1))", "sin(2*pi*x/12)"),
        x = c(x_ad, x_seas),
        value = round(c(r_ad["val"], r_seas["val"]), 6),
        ad_gradient = round(c(r_ad["grad"], r_seas["grad"]), 6),
        exact = round(c(expected_ad, expected_seas), 6),
        error = c(
            formatC(abs(r_ad["grad"] - expected_ad), format = "e", digits = 1),
            formatC(abs(r_seas["grad"] - expected_seas), format = "e", digits = 1)
        )
    ),
    row.names = FALSE,
    caption = "Dual number AD: automatic derivatives at machine precision"
)
```

## Caching

Compiled models are cached by a hash of the model specification string.
The cache directory is
`tools::R_user_dir("lucifer", "cache")/compiled_models/`. Recompilation
is skipped when a cached shared object exists for the same model hash,
reducing subsequent calls to
[`compile()`](https://robustecologies.github.io/lucifer/reference/compile.md)
to milliseconds. The cache is also used by parallel chain workers: each
child process loads the cached compiled object in ~100ms rather than
recompiling from source.

``` r

t1 <- system.time(compile(spec, Data, verbose = FALSE, cache = TRUE))[3]
t2 <- system.time(compile(spec, Data, verbose = FALSE, cache = TRUE))[3]

knitr::kable(
    data.frame(call = c("First", "Second (cached)"),
               seconds = round(c(t1, t2), 4)),
    row.names = FALSE,
    caption = "Compilation time with caching"
)
```

## Limitations

The C++ compilation backend has well-defined scope boundaries.

**Distribution coverage.** 24 of 86 registered distributions have C++
gradient support. Multivariate distributions (MVN, MVT, Wishart,
Dirichlet), specialized priors (HorseShoe, Lasso, Zellner), and exotic
distributions (Stable, VonMises, copulas) require the R backend. The
supported set covers standard GLMs, hierarchical normal models, and most
commonly used prior specifications.

**Chain rule depth.** The symbolic gradient generator handles one level
of function composition around matrix multiplication: `f(X %*% beta)`
where `f` is `exp`, `invlogit`, `log`, or `sqrt`. Deeper compositions
and conditional logic are not symbolically differentiated; the dual
number AD library provides the foundation for handling these cases but
is not yet integrated into the code generator.

**Posterior interpretation.** Since the compiled model operates in
unconstrained space, posterior samples require back-transformation:
[`exp()`](https://rdrr.io/r/base/Log.html) for positive parameters,
[`plogis()`](https://rdrr.io/r/stats/Logistic.html) for unit-interval
parameters. The sampler stores unconstrained values in the chain.

## Future directions

**AD integration into code generator.** The dual number library
(`lucifer/dual.h`) provides automatic differentiation for arbitrary
scalar expressions. Integrating it into the code generator would
eliminate the symbolic chain rule limitation: any deterministic
expression using standard arithmetic and elementary functions would be
differentiable automatically, including compositions like
`log(1 + exp(x))` (softplus), seasonal components (`sin(2*pi*t/12)`),
and custom link functions.

**Extended distribution coverage.** The gradient header already includes
`MVNGrad` and `DirichletGrad` structs. Wiring these into the code
generator requires extending the type dispatch to handle vector-valued
gradients and matrix-level chain rules, which would enable compiled
hierarchical multivariate models.

**Automatic posterior back-transformation.** Currently the user must
manually apply [`exp()`](https://rdrr.io/r/base/Log.html) or
[`plogis()`](https://rdrr.io/r/stats/Logistic.html) to posterior
samples. A future `as.constrained()` method could automate this based on
the `param_info` metadata stored in the compiled model.

## API reference

### `model_spec(formula, ..., backend)`

The `backend` argument accepts `"r"` (default), which generates an R
Model function compatible with all distributions and samplers; `"cpp"`,
which stores a flag so that
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
auto-compiles on first call when Data is available; and `"auto"`, which
attempts C++ compilation silently and falls back to R if compilation
fails.

### `compile(spec, Data, cache, verbose)`

Compiles a `model_spec` object to C++. Returns a `compiled_model` S3
object (inherits from `model_spec`) or the original spec on failure. Key
fields of the returned object include `$Model` (the C++ evaluation
function), `$Grad` (the analytical gradient function), `$initial_values`
(generating unconstrained-space starting values), `$cpp_code` (the full
generated source), `$hash` (8-character cache identifier), and `$n_parm`
(total parameter count).

### `code(compiled, section, file)`

The
[`code()`](https://robustecologies.github.io/lucifer/reference/code.md)
generic dispatches on `compiled_model` objects. The `section` argument
selects `"all"` (full source), `"eval"` (log-posterior function only),
`"grad"` (gradient function only), or `"data"` (data struct and
initialization). When `file` is provided, writes the source to disk and
returns the path invisibly.

### `check_nuts(fit, ebfmi_threshold, max_treedepth_frac, warmup)`

Runs NUTS-specific diagnostics on a `demonoid` fit object. Returns an S3
`nuts_check` object with [`print()`](https://rdrr.io/r/base/print.html)
and [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods.
Reports divergent transitions (count and percentage), E-BFMI per chain
(threshold default 0.3), and tree depth saturation (threshold default
0.1). For multi-chain fits, diagnostics are extracted from each chain’s
`$NUTS.Diagnostics` and overlaid in the plot with distinct colors. The
`warmup` argument (default 0) trims initial adaptation iterations from
diagnostic plots. The per-iteration data is available in
`fit$NUTS.Diagnostics` (single chain) or
`fit$Chains[[i]]$NUTS.Diagnostics` (multi-chain).

### Distribution gradient library

The `inst/include/lucifer/dist_gradients.h` header provides
`lucifer::grad::GradResult` structs with fields `lp` (log-pdf), `d_x`
(gradient w.r.t. random variable), `d_a1`, `d_a2`, `d_a3` (gradients
w.r.t. distribution parameters). These functions are available for
direct use in custom C++ code via `#include <lucifer/dist_gradients.h>`.

### Dual number AD library

The `inst/include/lucifer/dual.h` header provides `lucifer::ad::Dual`
structs for forward-mode automatic differentiation. A dual number
`Dual(value, derivative)` propagates derivatives through arithmetic and
20+ elementary functions. Include via `#include <lucifer/dual.h>` in
custom C++ code or sourceCpp scripts.

### Transform library

The `inst/include/lucifer/transforms.h` header provides
`log_transform()`, `logit_transform()`, `bounded_transform()`, and
`simplex_transform()` with corresponding inverse, log-Jacobian, and
gradient functions. These implement the same bijections used by Stan for
constrained parameter sampling.

## References

**\[1\]** Hoffman, M. D. and Gelman, A. (2014). The No-U-Turn Sampler:
adaptively setting path lengths in Hamiltonian Monte Carlo. *Journal of
Machine Learning Research*, 15, 1593–1623.

**\[2\]** Carpenter, B., Gelman, A., Hoffman, M. D., Lee, D., Goodrich,
B., Betancourt, M., Brubaker, M., Guo, J., Li, P. and Riddell, A.
(2017). Stan: a probabilistic programming language. *Journal of
Statistical Software*, 76(1).
[doi:10.18637/jss.v076.i01](https://doi.org/10.18637/jss.v076.i01).

**\[3\]** Eddelbuettel, D. and Francois, R. (2011). Rcpp: seamless R and
C++ integration. *Journal of Statistical Software*, 40(8), 1–18.
[doi:10.18637/jss.v040.i08](https://doi.org/10.18637/jss.v040.i08).

**\[4\]** de Valpine, P., Turek, D., Paciorek, C. J., Anderson-Bergman,
C., Lang, D. T. and Bodik, R. (2017). Programming with models: writing
statistical algorithms for general model structures with NIMBLE.
*Journal of Computational and Graphical Statistics*, 26(2), 403–413.
[doi:10.1080/10618600.2016.1172487](https://doi.org/10.1080/10618600.2016.1172487).

**\[5\]** Kristensen, K., Nielsen, A., Berg, C. W., Skaug, H. and Bell,
B. M. (2016). TMB: automatic differentiation and Laplace approximation.
*Journal of Statistical Software*, 70(5), 1–21.
[doi:10.18637/jss.v070.i05](https://doi.org/10.18637/jss.v070.i05).

**\[6\]** Betancourt, M. (2017). A conceptual introduction to
Hamiltonian Monte Carlo. *arXiv:1701.02434*.

**\[7\]** Betancourt, M. (2016). Diagnosing suboptimal cotangent
disintegrations in Hamiltonian Monte Carlo. *arXiv:1604.00695*.
