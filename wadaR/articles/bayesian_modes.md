# Bayesian mode inference with lucifer

## Why modes matter more than means

Bayesian inference typically summarizes posterior distributions through
their mean or median, but this practice can be dangerously misleading
when posteriors are multimodal. Consider a bimodal distribution with
peaks at 0 and 10: the posterior mean of 5 falls in a region of
near-zero density, representing a value that the data never actually
support. In decision-theoretic terms, the mode (maximum a posteriori,
MAP) minimizes 0-1 loss, the mean minimizes squared-error loss, and the
median minimizes absolute-error loss. When the goal is to identify the
single most probable value, the mode is the only appropriate point
estimate.

Multimodal posteriors arise more frequently than practitioners expect.
Mixture models, hierarchical specifications with weak identifiability,
and models with discrete latent structures routinely produce posterior
surfaces with multiple peaks. Detecting and characterizing these modes
is essential for understanding model behavior, yet most software reports
only means and credible intervals, silently averaging over qualitatively
different posterior regions.

The problem extends beyond Bayesian inference. Any empirical
distribution, whether a posterior sample, a dataset of economic
observations, or a vector of biological measurements, may contain
multiple modes that reveal distinct subpopulations or regimes. Classical
tests for multimodality (Silverman’s bandwidth test, Hartigan’s dip
test) are frequentist and provide only binary reject/fail-to-reject
decisions. Bayesian mode inference, by contrast, produces a full
posterior distribution over both the number and location of modes,
enabling nuanced probabilistic statements.

``` r

x <- c(rnorm(500, 0, 1), rnorm(500, 10, 1))
m <- Modes(x)
dat <- data.frame(x = x)
mode_df <- data.frame(xint = m$modes, label = "Mode")
mean_df <- data.frame(xint = mean(x), label = "Mean")
ggplot2::ggplot(dat, ggplot2::aes(x = x)) +
  ggplot2::geom_density(fill = "#4e79a7", alpha = 0.3, color = "#4e79a7",
                        linewidth = 0.8) +
  ggplot2::geom_vline(data = mean_df, ggplot2::aes(xintercept = xint),
                      color = "#e15759", linewidth = 1, linetype = "dashed") +
  ggplot2::geom_vline(data = mode_df, ggplot2::aes(xintercept = xint),
                      color = "#59a14f", linewidth = 1) +
  ggplot2::annotate("text", x = mean(x), y = Inf, label = "Mean",
                    vjust = 2, color = "#e15759", fontface = "bold", size = 3.5) +
  ggplot2::annotate("text", x = m$modes, y = Inf,
                    label = paste0("Mode ", seq_along(m$modes)),
                    vjust = 3.5, color = "#59a14f", fontface = "bold", size = 3.5) +
  theme_relab() +
  ggplot2::labs(x = "Value", y = "Density",
                title = "Bimodal posterior: mean vs. modes",
                subtitle = "The mean (dashed red) falls in a low-density region")
```

## Classical mode estimation

lucifer provides two primary mode estimation functions. `Mode(x)`
returns the single most prominent mode using kernel density estimation
for continuous data or frequency counting for discrete data. `Modes(x)`
returns all detected modes ordered by density, each with its relative
size (the proportion of the density it captures).

Both functions support four estimation methods: kernel density
estimation (`"kde"`, the default), the half-sample mode (`"hsm"`;
[\[4\]](#ref4)), Venter’s mode (`"venter"`), and the shortest half-range
(`"shorth"`). The half-sample mode is particularly robust to outliers
since it recursively selects the densest half of the sample without
requiring bandwidth selection. Venter’s mode finds the shortest interval
containing 50% of the data and returns its midpoint, providing a
different notion of centrality.

``` r

set.seed(42)
x <- c(rnorm(400, -2, 0.8), rnorm(300, 3, 1.2), rnorm(300, 8, 0.6))
ms <- Modes(x)
print(ms)
```

``` r

summary(ms, compare = TRUE)
```

``` r

plot(ms, type = "density")
```

The modality testing functions
[`is.unimodal()`](https://robustecologies.github.io/lucifer/reference/Mode.md),
[`is.bimodal()`](https://robustecologies.github.io/lucifer/reference/Mode.md),
[`is.trimodal()`](https://robustecologies.github.io/lucifer/reference/Mode.md),
and
[`is.multimodal()`](https://robustecologies.github.io/lucifer/reference/Mode.md)
provide quick Boolean assessments:

``` r

is.unimodal(rnorm(1000))
is.trimodal(x)
```

## The Bayesian approach: Sparse Finite Mixtures

Classical mode detection treats the distribution as fixed and known (up
to a bandwidth parameter). The Bayesian approach, developed by Cross,
Hoogerheide, Ulker and van Dijk [\[1\]](#ref1), instead fits a mixture
model to the data and then extracts modes from the fitted mixture at
each MCMC iteration, producing a posterior distribution over the modes
themselves. This yields three quantities that classical methods cannot
provide: the posterior probability that the distribution is multimodal,
a posterior distribution over the number of modes, and credible
intervals for each mode location.

### The Sparse Finite Mixture framework

The method uses Sparse Finite Mixtures (SFM; [\[2\]](#ref2)) to
simultaneously fit the mixture and estimate the number of components.
The key idea is to deliberately overspecify the number of components
\\K\\ and let a sparsity-inducing prior on the Dirichlet concentration
parameter \\\alpha\\ shrink unnecessary components to zero weight.

The mixture model for data \\y_1, \ldots, y_n\\ is

\\p(y_i \mid \boldsymbol{\theta}, K) = \sum\_{k=1}^{K} \pi_k \\ f(y_i
\mid \theta_k)\\

where \\\pi_k\\ are the mixture weights and \\f(\cdot \mid \theta_k)\\
is the component density with parameters \\\theta_k\\. The weights
follow a symmetric Dirichlet prior \\(\pi_1, \ldots, \pi_K) \sim
\text{Dirichlet}(\alpha, \ldots, \alpha)\\, and the concentration
parameter has a Gamma hyperprior \\\alpha \sim \text{Gamma}(a\_\alpha,
b\_\alpha)\\ with \\E(\alpha) = a\_\alpha / b\_\alpha\\ set to a small
value (default \\1/200\\) that strongly favors parsimony.

### Two-stage algorithm

**Stage 1 (mixture fitting):** A Gibbs MCMC sampler iterates over:

1.  Component allocations \\z_i \sim \text{Multinomial}(\pi_1 f(y_i \mid
    \theta_1), \ldots, \pi_K f(y_i \mid \theta_K))\\
2.  Mixture weights \\\boldsymbol{\pi} \sim \text{Dirichlet}(\alpha +
    n_1, \ldots, \alpha + n_K)\\
3.  Concentration \\\alpha\\ via Metropolis-Hastings
4.  Component parameters \\\theta_k\\ from their family-specific
    posteriors

**Stage 2 (mode detection):** At each post-burn-in MCMC draw, modes of
the fitted mixture density are detected using:

- **Fixed-point iteration** for Normal mixtures ([\[3\]](#ref3)):
  converges to modes from each component mean
- **Discrete enumeration** for count distributions: evaluates the
  mixture PMF on integer support and finds peaks
- **Grid search with golden-section refinement** for other continuous
  distributions

The result is a collection of mode sets, one per MCMC draw, from which
the posterior distribution over mode count and locations is computed
directly.

## Available distribution families

lucifer implements 11 mixture component families for
[`BayesMode()`](https://robustecologies.github.io/lucifer/reference/BayesMode.md):

When `family = "auto"` (the default), lucifer selects based on data
properties: integer non-negative data maps to `"poisson"` (or
`"shifted_poisson"` if all values are positive), positive continuous
data to `"gamma"`, and general continuous data to `"normal"`.

## Tutorial: continuous bimodal data

``` r

set.seed(42)
x <- c(rnorm(400, 0, 1), rnorm(400, 7, 1.5))
bm <- BayesMode(x, family = "normal", K = 10, iter = 5000, burnin = 2000)
```

The print method provides a clean summary:

``` r

print(bm)
```

The summary includes MCMC diagnostics:

``` r

summary(bm)
```

### Visualizations

The default plot shows the data density with posterior mode locations as
colored vertical bands whose width represents the 95% credible interval
and opacity reflects the posterior probability:

``` r

plot(bm, type = "modes")
```

The number-of-modes posterior:

``` r

plot(bm, type = "n_modes")
```

The fitted mixture density overlaid on the data, showing both individual
posterior draws and the posterior mean:

``` r

plot(bm, type = "mixture")
```

MCMC trace plots for diagnostic quantities:

``` r

plot(bm, type = "trace")
```

Violin plots showing the full posterior distribution of mode locations:

``` r

plot(bm, type = "posterior")
```

Convergence diagnostics via running means:

``` r

plot(bm, type = "convergence")
```

## Tutorial: discrete count data

Bayesian mode inference is particularly valuable for count data, where
classical kernel-based methods are inappropriate. The shifted-Poisson
family from Cross et al. [\[1\]](#ref1) handles both under- and
overdispersion through the shift parameter \\\kappa_k\\, which decouples
the mean from the variance within each component.

``` r

set.seed(42)
x <- c(rpois(300, 3), rpois(300, 12))
bm_pois <- BayesMode(x, family = "poisson", K = 8,
                      iter = 5000, burnin = 2000)
print(bm_pois)
```

``` r

plot(bm_pois, type = "modes")
```

``` r

plot(bm_pois, type = "n_modes")
```

## Tutorial: posterior parameter modes

When passed a fitted lucifer object (class `demonoid`, `laplace`,
`iterquad`, `vb`, or `pmc`),
[`BayesMode()`](https://robustecologies.github.io/lucifer/reference/BayesMode.md)
automatically extracts the posterior samples for each parameter and
performs mode inference independently on each marginal posterior. This
is invaluable for detecting multimodal posteriors that may indicate
model identifiability issues or genuine multi-modality in the parameter
space.

``` r

# After fitting a model with lucifer():
# result <- lucifer(Model, Data, Initial.Values,
#                   Covar = NULL, Iterations = 10000,
#                   Algorithm = "NUTS")
# bm_post <- BayesMode(result)
# print(bm_post)
# plot(bm_post)
```

The `bayes_mode_multi` print method provides a compact
one-line-per-parameter summary showing the MAP number of modes,
P(unimodal), and the primary mode location with its credible interval.

## Advanced: tuning and diagnostics

### Choosing K

The SFM framework is robust to the choice of \\K\\ because the sparsity
prior actively prunes unnecessary components. In practice, \\K = 10\\
works well for most datasets; increase to \\K = 15\\–\\20\\ if you
suspect more than 5 modes. Setting \\K\\ too low risks missing modes;
setting it too high increases computation but does not bias the results.

### The alpha prior

The Gamma hyperprior on \\\alpha\\ controls the sparsity. The default
\\\text{Gamma}(1, 200)\\ gives \\E(\alpha) = 0.005\\, strongly favoring
few active components. For exploratory analysis where you want to err on
the side of detecting modes, use `alpha.prior = c(1, 50)` for a less
sparse prior. For confirmatory analysis where parsimony matters, the
default is appropriate.

### MCMC convergence

The summary method reports ESS for the Dirichlet concentration
\\\alpha\\ and the mode count. Low ESS for \\\alpha\\ (below 50)
suggests poor mixing of the concentration parameter; increase `iter` or
adjust the MH proposal. The trace plots (`plot(bm, type = "trace")`) and
convergence plots (`plot(bm, type = "convergence")`) provide visual
diagnostics.

### Computational performance

The C++ engine with OpenMP parallelization handles datasets of several
thousand observations efficiently. For \\n = 1000\\, \\K = 10\\, and
5000 iterations, typical runtime is 1–3 seconds on a modern machine. The
allocation step (which scales as \\O(n \times K)\\) is the computational
bottleneck.

## Mathematical appendix

### Full conditionals for Normal mixture

For a Normal mixture with \\f(y \mid \mu_k, \sigma_k) = \mathcal{N}(y;
\mu_k, \sigma_k^2)\\:

**Weights:** \\\boldsymbol{\pi} \mid \mathbf{z} \sim
\text{Dirichlet}(\alpha + n_1, \ldots, \alpha + n_K)\\

**Mean:** \\\mu_k \mid \sigma_k, \mathbf{y}, \mathbf{z} \sim
\mathcal{N}\left(\frac{b_0 / B_0 + \bar{y}\_k n_k / \sigma_k^2}{1/B_0 +
n_k / \sigma_k^2}, \\ \frac{1}{1/B_0 + n_k/\sigma_k^2}\right)\\

**Variance:** \\\sigma_k^{-2} \mid \mu_k, \mathbf{y}, \mathbf{z} \sim
\text{Gamma}\left(c_0 + \frac{n_k}{2}, \\ C_0 + \frac{1}{2}\sum\_{i: z_i
= k}(y_i - \mu_k)^2\right)\\

### Full conditionals for Poisson mixture

For a Poisson mixture with \\f(y \mid \lambda_k) = \text{Pois}(y;
\lambda_k)\\:

**Rate:** \\\lambda_k \mid \mathbf{y}, \mathbf{z} \sim
\text{Gamma}\left(l_0 + \sum\_{i: z_i=k} y_i, \\ L_0 + n_k\right)\\

### Mode detection: fixed-point iteration

For a Normal mixture, modes satisfy \\\nabla \log f(x) = 0\\. The
fixed-point iteration

\\x^{(t+1)} = \frac{\sum\_{k=1}^K p(k \mid x^{(t)}) \mu_k /
\sigma_k^2}{\sum\_{k=1}^K p(k \mid x^{(t)}) / \sigma_k^2}\\

where \\p(k \mid x) = \pi_k \phi(x; \mu_k, \sigma_k) / \sum_j \pi_j
\phi(x; \mu_j, \sigma_j)\\, converges to a mode when initialized from
each component mean.

### Mode definition for discrete distributions

A mode \\y_m\\ of a discrete mixture satisfies \\f(y_m - 1) \< f(y_m) \>
f(y_m + 1)\\ for interior modes, or \\f(0) \> f(1)\\ for a boundary mode
at zero. Plateau modes (consecutive equal-probability values) are
counted as a single mode.

## References

**\[1\]** Cross, J.L., Hoogerheide, L., Ulker, S. and van Dijk, H.K.
(2024). Bayesian mode inference for discrete distributions in economics
and finance. *Economics Letters*, 235, 111579.
[doi:10.1016/j.econlet.2024.111579](https://doi.org/10.1016/j.econlet.2024.111579)

**\[2\]** Malsiner-Walli, G., Fruhwirth-Schnatter, S. and Grun, B.
(2016). Model-based clustering based on sparse finite Gaussian mixtures.
*Statistics and Computing*, 26, 303-324.
[doi:10.1007/s11222-014-9500-2](https://doi.org/10.1007/s11222-014-9500-2)

**\[3\]** Carreira-Perpinán, M.A. (2000). Mode-finding for mixtures of
Gaussian distributions. *IEEE Transactions on Pattern Analysis and
Machine Intelligence*, 22(11), 1318-1323.
[doi:10.1109/34.888718](https://doi.org/10.1109/34.888718)

**\[4\]** Bickel, D.R. and Fruhwirth, R. (2006). On a fast, robust
estimator of the mode: comparisons to other robust estimators with
applications. *Computational Statistics & Data Analysis*, 50(12),
3500-3530.
[doi:10.1016/j.csda.2005.07.014](https://doi.org/10.1016/j.csda.2005.07.014)
