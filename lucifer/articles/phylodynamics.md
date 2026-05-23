# Phylodynamic inference with lucifer

## Introduction

Bayesian phylogenetic and phylodynamic inference combines molecular
sequence data with evolutionary models to reconstruct the history of
lineages through time. The dominant software for this task, BEAST
[\[1\]](#ref1) [\[2\]](#ref2), jointly samples tree topologies and
continuous parameters (substitution rates, clock rates, coalescent
demographic parameters) using a specialized MCMC engine with
tree-specific proposal mechanisms such as subtree-pruning-and-regrafting
(SPR) and nearest-neighbor interchange (NNI) moves.

lucifer is not a tree sampler. It has no SPR proposals, no NNI moves, no
machinery for traversing discrete tree space. What it does have is a
flexible Bayesian inference engine with over 40 MCMC algorithms,
sequential Monte Carlo (SMC) with marginal likelihood estimation, modern
convergence diagnostics (rank-normalized \\\hat{R}\\, bulk and tail
ESS), leave-one-out cross-validation via Pareto smoothed importance
sampling (LOO-PSIS), Bayesian model stacking, and automated algorithm
selection through the Crucible pipeline. These capabilities offer
genuine complementary value for phylodynamic problems when the tree
topology is fixed.

Conditioning on a fixed tree is not a toy exercise. In practice,
researchers routinely fix trees obtained from maximum likelihood (e.g.,
RAxML, IQ-TREE) or consensus methods and estimate downstream parameters.
Gao et al. [\[3\]](#ref3) recently demonstrated that tree-topology
mixing is precisely where standard BEAST diagnostics fail: posterior
tree landscapes are inherently multimodal, with multiple peaks separated
by deep valleys that MCMC chains struggle to cross even after billions
of iterations. Their analysis of 15 large phylodynamic datasets showed
that standard continuous-parameter diagnostics (ESS, \\\hat{R}\\) often
fail to detect these tree-space sampling problems, which in turn distort
estimates of evolutionary rates, divergence times, and demographic
trajectories. In this context, explicitly conditioning on tree
hypotheses and comparing parameter estimates across them is arguably
more transparent than relying on a joint sampler whose topology mixing
may be inadequate.

This vignette demonstrates how lucifer addresses core phylodynamic
inference problems using its existing infrastructure. We cover five
examples: (1) GTR+\\\Gamma\\ substitution model estimation on a fixed
tree; (2) model comparison among competing substitution models via SMC
marginal likelihoods and LOO-PSIS; (3) strict versus relaxed molecular
clock estimation as a hierarchical model; (4) coalescent demographic
inference including Bayesian skyline reconstruction; and (5)
multimodality detection and diagnostics connecting to the findings of
Gao et al. [\[3\]](#ref3). All examples use data from the `ape`
[\[4\]](#ref4) and `phangorn` [\[5\]](#ref5) R packages, requiring no
external downloads.

``` r

library(lucifer)
library(ape)
library(phangorn)
library(ggplot2)
```

## Background theory

### Substitution models

The evolution of nucleotide sequences along a phylogenetic tree is
modeled as a continuous-time Markov chain (CTMC) on the state space
\\\\A, C, G, T\\\\. The instantaneous rate matrix \\Q\\ governs the
transition rates between nucleotides; its off-diagonal entries
\\q\_{ij}\\ represent the instantaneous rate of substitution from
nucleotide \\i\\ to \\j\\, and its diagonal entries satisfy \\q\_{ii} =
-\sum\_{j \neq i} q\_{ij}\\ so that each row sums to zero. Given a
branch of length \\t\\ (measured in expected substitutions per site),
the transition probability matrix is the matrix exponential \\P(t) =
\exp(Qt)\\, where entry \\P\_{ij}(t)\\ gives the probability that
nucleotide \\i\\ at the ancestral end of the branch becomes nucleotide
\\j\\ at the descendant end.

The general time-reversible (GTR) model [\[6\]](#ref6) parameterizes
\\Q\\ as \\Q\_{ij} = \rho\_{ij} \pi_j\\ for \\i \neq j\\, where
\\\boldsymbol{\pi} = (\pi_A, \pi_C, \pi_G, \pi_T)\\ are the stationary
(equilibrium) base frequencies satisfying \\\sum_k \pi_k = 1\\, and
\\\rho\_{ij} = \rho\_{ji}\\ are the symmetric exchangeability
parameters. The six upper-triangular exchangeabilities are
conventionally written as \\(a, b, c, d, e, f)\\ corresponding to the
pairs (AC, AG, AT, CG, CT, GT), with \\f = 1\\ fixed for
identifiability. The GTR model thus has 8 free parameters: 5
exchangeability rates and 3 free base frequencies (the fourth being
determined by the simplex constraint).

Two important special cases reduce the parameter count. The Jukes-Cantor
(JC69) model [\[7\]](#ref7) sets all exchangeabilities equal and all
base frequencies to 0.25, yielding zero free substitution parameters.
The Hasegawa-Kishino-Yano (HKY) model [\[8\]](#ref8) distinguishes
transitions (purine\\\leftrightarrow\\purine and
pyrimidine\\\leftrightarrow\\pyrimidine, governed by \\\kappa\\) from
transversions (all other substitutions, rate 1), with empirical base
frequencies, giving one free substitution parameter.

Rate heterogeneity across sites is modeled by a discrete gamma
distribution [\[9\]](#ref9) with shape parameter \\\alpha\\ and \\K\\
rate categories (typically \\K = 4\\). When \\\alpha\\ is small, rates
vary substantially across sites; as \\\alpha \to \infty\\, all sites
evolve at the same rate. This adds one free parameter to any
substitution model.

The likelihood of an alignment \\D\\ given a tree \\T\\ and model
parameters \\\theta\\ is computed by Felsenstein’s pruning algorithm
[\[10\]](#ref10), which traverses the tree from tips to root, combining
partial likelihood vectors at each internal node. For \\N\\ alignment
sites with \\K\\ unique patterns and pattern weights \\w_k\\, the total
log-likelihood decomposes as \\\log L(D \mid T, \theta) =
\sum\_{k=1}^{K} w_k \log L_k\\, where \\L_k\\ is the likelihood of
pattern \\k\\. This per-pattern decomposition is critical for LOO-PSIS,
which requires per-observation log-likelihoods.

### Molecular clock models

A molecular clock relates the substitution branch lengths (in expected
substitutions per site) to divergence times (in calendar units) through
evolutionary rates. Under the strict clock [\[11\]](#ref11), a single
global rate \\r\\ applies to all branches, so the substitution branch
length for branch \\i\\ is \\b_i = r \cdot t_i\\, where \\t_i\\ is the
duration of branch \\i\\ in time units. This assumption is restrictive:
empirical data frequently show rate variation across lineages.

The uncorrelated relaxed lognormal clock [\[12\]](#ref12) assigns each
branch its own rate \\r_i\\ drawn independently from a lognormal
distribution: \\r_i \sim \text{LogNormal}(\mu_r, \sigma_r)\\. The
substitution branch length becomes \\b_i = r_i \cdot t_i\\. This is a
hierarchical model: the branch-specific rates are partially pooled
toward the population mean through the shared hyperparameters \\(\mu_r,
\sigma_r)\\. When \\\sigma_r \to 0\\, the relaxed clock collapses to a
strict clock; when \\\sigma_r\\ is large, branches evolve at very
different rates. This hierarchical structure maps naturally onto
lucifer’s `Model()` interface, where the branch rates are parameters
with lognormal priors governed by estimated hyperparameters.

### Coalescent demography

The coalescent [\[13\]](#ref13) provides a probabilistic framework for
relating a genealogy (the tree of sampled lineages) to the demographic
history of the population from which they were sampled. Given \\n\\
sampled lineages, the process works backward in time: at any point where
\\k\\ lineages coexist, the waiting time until the next coalescence
event is exponentially distributed with rate \\\binom{k}{2} / N_e(t)\\,
where \\N_e(t)\\ is the effective population size at time \\t\\.

For a rooted ultrametric tree with \\n\\ tips and \\n - 1\\ coalescence
events at times \\t_1 \< t_2 \< \ldots \< t\_{n-1}\\ (measured backward
from the present), the coalescent log-likelihood given a demographic
function \\N_e(t)\\ is

\\\log L(T \mid N_e) = \sum\_{j=1}^{n-1} \left\[ \log \binom{k_j}{2} -
\log N_e(t_j) - \int\_{t\_{j-1}}^{t_j} \frac{\binom{k_j}{2}}{N_e(s)} \\
ds \right\]\\

where \\k_j\\ is the number of lineages in interval \\j\\ (between
coalescence events \\j-1\\ and \\j\\, with \\t_0 = 0\\), and the
integral represents the cumulative hazard over that interval.

Three demographic parameterizations are commonly used. Under **constant
population size**, \\N_e(t) = N_e\\ for all \\t\\, and the integral
simplifies to \\\binom{k_j}{2} \Delta t_j / N_e\\ where \\\Delta t_j =
t_j - t\_{j-1}\\. Under **exponential growth**, \\N_e(t) = N_0
\exp(gt)\\ where \\g\\ is the growth rate and \\N_0\\ is the present-day
size; the integral has a closed-form solution involving the exponential
function. The **Bayesian skyline** [\[14\]](#ref14) models \\N_e\\ as
piecewise constant over grouped coalescent intervals, with \\M\\
population size parameters \\N_1, \ldots, N_M\\, each governing a
contiguous block of coalescent intervals. This non-parametric approach
can recover complex demographic trajectories without imposing a
functional form.

All three parameterizations yield analytical coalescent log-likelihoods
that can be coded directly in R, requiring no external phylogenetic
library.

## Preparing the data

We use the `woodmouse` dataset from the `ape` package: 15 European
woodmouse (*Apodemus sylvaticus*) mitochondrial cytochrome *b* gene
sequences, 965 base pairs, published by Michaux et al. [\[15\]](#ref15).
We construct a neighbor-joining tree and root it for downstream
analyses.

``` r

# Load sequence data
data(woodmouse, package = "ape")

# Compute distance matrix and build NJ tree
dm <- ape::dist.dna(woodmouse, model = "F84")
tree_nj <- ape::nj(dm)

# Root the tree (required for clock models and coalescent)
tree_rooted <- ape::root(tree_nj, outgroup = "No305", resolve.root = TRUE)

# Convert sequences to phangorn format
pdat <- phangorn::as.phyDat(woodmouse)

cat("Sequences:", length(labels(woodmouse)), "\n")
cat("Sites:", ncol(woodmouse), "\n")
cat("Unique patterns:", attr(pdat, "nr"), "\n")
cat("Tree edges:", nrow(tree_rooted$edge), "\n")
```

### The phylogenetic likelihood helper

We define a small helper function that wraps
[`phangorn::pml()`](https://klausvigo.github.io/phangorn/reference/pml.html)
to accept a lucifer-style parameter vector. This function extracts the
base frequencies, GTR exchangeability rates, and gamma shape parameter
from the parameter vector, calls
[`pml()`](https://klausvigo.github.io/phangorn/reference/pml.html), and
returns the log-likelihood along with per-site log-likelihoods for
LOO-PSIS. Since
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
runs parallel chains in separate R processes via `callr`, helper
functions defined in the calling environment are not accessible from
within the `Model()` function. The standard workaround is to pass the
helper through the `Data` list, which is serialized and transmitted to
each subprocess. All `Data` lists in this vignette therefore include a
`phylo_loglik` element.

``` r

#' Compute phylogenetic log-likelihood via Felsenstein pruning
#' @param bf Numeric vector of 4 base frequencies (must sum to 1)
#' @param Q Numeric vector of 6 GTR exchangeability rates
#' @param shape Gamma shape parameter for rate heterogeneity
#' @param tree A phylo object (the fixed tree)
#' @param pdat A phyDat object (the sequence alignment)
#' @param k Number of gamma rate categories (default 4)
#' @return List with logLik (scalar) and site_ll (per-site log-likelihoods)
#'
#' The closure includes a gc() call every 100 pml() evaluations. phangorn's
#' C-level likelihood engine accumulates internal allocations that can
#' corrupt the heap when called hundreds of thousands of times in a single
#' R process; periodic garbage collection keeps the allocations compact.
phylo_loglik <- local({
    n <- 0L
    function(bf, Q, shape, tree, pdat, k = 4) {
        n <<- n + 1L
        if (n %% 100L == 0L) gc(verbose = FALSE)
        fit <- phangorn::pml(tree, pdat, bf = bf, Q = Q, k = k, shape = shape)
        site_ll <- rep(fit$siteLik, attr(pdat, "weight"))
        list(logLik = fit$logLik, site_ll = site_ll)
    }
})
```

## Example 1: GTR+\\\Gamma\\ substitution model estimation

This example estimates the full GTR+\\\Gamma\\ substitution model
parameters on the fixed NJ tree. The parameter vector has 9 components:
3 free base frequencies (the fourth determined by the simplex
constraint), 5 free GTR exchangeability rates (the sixth, GT, fixed at 1
for identifiability), and 1 gamma shape parameter. Base frequencies are
parameterized through a softmax (additive log-ratio) transformation: the
sampler operates on three unconstrained real-valued parameters
\\\alpha_A, \alpha_C, \alpha_G\\ (with \\\alpha_T = 0\\ as reference),
and the frequencies are recovered as \\\pi_k = \exp(\alpha_k) / \sum_j
\exp(\alpha_j)\\. This avoids the cascading
[`interval()`](https://robustecologies.github.io/lucifer/reference/interval.md)
constraints that would otherwise break componentwise samplers such as
AMWG when one frequency is proposed near its boundary.

### Model specification

``` r

# ---- Data list ----
n_sites <- ncol(woodmouse)

Data <- list(
    # Required lucifer fields
    mon.names  = c("LP", "LL"),
    parm.names = c("alpha_A", "alpha_C", "alpha_G",
                   "q_AC", "q_AG", "q_AT", "q_CG", "q_CT",
                   "shape"),
    # Phylogenetic data (passed through to Model)
    tree = tree_rooted,
    pdat = pdat,
    N    = n_sites,
    phylo_loglik = phylo_loglik
)

# ---- Model function ----
Model <- function(parm, Data) {

    # --- Base frequencies via softmax (unconstrained parameterization) ---
    # parm[1:3] are log-ratio parameters (alpha_A, alpha_C, alpha_G);
    # alpha_T = 0 is the reference category. This avoids cascading
    # interval() constraints that break componentwise samplers like AMWG.
    alpha <- c(parm[1:3], 0)
    alpha <- alpha - max(alpha)  # numerical stability
    bf    <- exp(alpha) / sum(exp(alpha))

    # GTR rates: all positive, GT fixed at 1
    q_AC <- interval(parm[4], 1e-4, 100)
    q_AG <- interval(parm[5], 1e-4, 100)
    q_AT <- interval(parm[6], 1e-4, 100)
    q_CG <- interval(parm[7], 1e-4, 100)
    q_CT <- interval(parm[8], 1e-4, 100)
    parm[4:8] <- c(q_AC, q_AG, q_AT, q_CG, q_CT)

    Q <- c(q_AC, q_AG, q_AT, q_CG, q_CT, 1)  # GT = 1 (reference)

    # Gamma shape
    shape <- interval(parm[9], 0.01, 100)
    parm[9] <- shape

    # --- Log-likelihood ---
    res <- Data$phylo_loglik(bf, Q, shape, Data$tree, Data$pdat)
    LL  <- res$logLik

    # --- Priors ---
    # Dirichlet(1,1,1,1) is uniform on the simplex; in the unconstrained
    # alpha space the prior contributes only the log-Jacobian of the
    # softmax transformation: log|J| = sum(log(pi_k))
    pi_prior <- sum(log(bf))
    # Exponential(1) on each GTR rate
    q_prior <- sum(dexp(c(q_AC, q_AG, q_AT, q_CG, q_CT),
                        rate = 1, log = TRUE))
    # Exponential(1) on gamma shape
    shape_prior <- dexp(shape, rate = 1, log = TRUE)

    LP <- LL + pi_prior + q_prior + shape_prior

    # --- Posterior predictive (monitor log-likelihood) ---
    yhat <- rep(LL / Data$N, Data$N)

    list(LP = LP, Dev = -2 * LL, Monitor = c(LP, LL),
         yhat = yhat, parm = parm)
}

# ---- Initial values ----
# Start near uniform/default values
# First 3 are softmax log-ratios: (0,0,0) maps to uniform (0.25 each)
Initial.Values <- c(
    0, 0, 0,                 # softmax alphas (uniform start)
    1, 1, 1, 1, 1,          # GTR rates (all equal = JC-like start)
    1                        # gamma shape
)
```

### Fitting with AMWG

The phylogenetic likelihood is computed through
[`phangorn::pml()`](https://klausvigo.github.io/phangorn/reference/pml.html),
an opaque compiled function for which analytical gradients are not
available. Gradient-based samplers such as NUTS would resort to
finite-difference approximations, requiring \\2p\\ additional likelihood
evaluations per gradient (18 for this 9-parameter model) and many
gradients per leapfrog trajectory, making them prohibitively slow.
Prescribe confirms this diagnosis: it detects that gradients are
unavailable and recommends componentwise Metropolis algorithms. We use
AMWG (Adaptive Metropolis-within-Gibbs [\[20\]](#ref20)), which updates
each parameter individually with automatically tuned proposal scales and
requires only one likelihood evaluation per component update.

``` r

# Run Prescribe to verify algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

# Fit with AMWG (componentwise; no gradients needed)
fit_gtr <- lucifer(
    Model, Data, Initial.Values,
    Iterations = 20000, Status = 5000, Thinning = 10,
    Algorithm = "AMWG",
    Specs = list(B = NULL, n = 0, Periodicity = 50)
)

# Diagnostics
print(fit_gtr)
Consort(fit_gtr)
```

### Posterior summaries

``` r

# Extract posterior samples (after burnin)
post <- fit_gtr$Posterior2
cat("Posterior dimensions:", dim(post), "\n")

# Transform softmax alphas to base frequencies for reporting
bf_post <- t(apply(post[, 1:3, drop = FALSE], 1, function(a) {
    alpha <- c(a, 0)
    alpha <- alpha - max(alpha)
    exp(alpha) / sum(exp(alpha))
}))
colnames(bf_post) <- c("pi_A", "pi_C", "pi_G", "pi_T")

cat("\n--- Base frequency posteriors ---\n")
for (i in 1:4) {
    cat(sprintf("  %s: mean=%.4f, sd=%.4f, 95%% CI=[%.4f, %.4f]\n",
                colnames(bf_post)[i], mean(bf_post[, i]), sd(bf_post[, i]),
                quantile(bf_post[, i], 0.025), quantile(bf_post[, i], 0.975)))
}

cat("\n--- GTR rate posteriors ---\n")
for (i in 4:8) {
    cat(sprintf("  %s: mean=%.4f, sd=%.4f, 95%% CI=[%.4f, %.4f]\n",
                Data$parm.names[i], mean(post[, i]), sd(post[, i]),
                quantile(post[, i], 0.025), quantile(post[, i], 0.975)))
}

cat("\n--- Gamma shape posterior ---\n")
cat(sprintf("  shape: mean=%.4f, sd=%.4f, 95%% CI=[%.4f, %.4f]\n",
            mean(post[, 9]), sd(post[, 9]),
            quantile(post[, 9], 0.025), quantile(post[, 9], 0.975)))

# Compare with ML estimates
fit_ml <- phangorn::optim.pml(
    phangorn::pml(tree_rooted, pdat), model = "GTR",
    optGamma = TRUE, rearrangement = "none",
    control = phangorn::pml.control(trace = 0))

cat("\n--- ML comparison ---\n")
cat("ML bf:", round(fit_ml$bf, 4), "\n")
cat("ML Q:", round(fit_ml$Q, 4), "\n")
cat("ML shape:", round(fit_ml$shape, 4), "\n")
cat("ML logLik:", round(fit_ml$logLik, 2), "\n")
```

### Convergence diagnostics

``` r

# Modern diagnostics
rhat_vals <- sapply(1:9, function(i) Rhat(post[, i]))
ess_bulk  <- sapply(1:9, function(i) ESS.bulk(post[, i]))
ess_tail  <- sapply(1:9, function(i) ESS.tail(post[, i]))

diag_df <- data.frame(
    Parameter = Data$parm.names,
    Rhat      = round(rhat_vals, 4),
    ESS.bulk  = round(ess_bulk, 0),
    ESS.tail  = round(ess_tail, 0)
)
knitr::kable(diag_df, caption = "Convergence diagnostics for GTR+G model")
```

## Example 2: model comparison among substitution models

A fundamental question in phylogenetics is which substitution model best
describes the data. We compare three models of increasing complexity on
the same woodmouse alignment and fixed tree: JC69 (0 free substitution
parameters), HKY (1 free parameter: \\\kappa\\), and GTR (8 free
parameters). Each model includes gamma rate heterogeneity with 4
categories.

### JC69+\\\Gamma\\ model

``` r

Data_jc <- list(
    mon.names  = c("LP", "LL"),
    parm.names = c("shape"),
    tree = tree_rooted,
    pdat = pdat,
    N    = n_sites,
    phylo_loglik = phylo_loglik
)

Model_jc <- function(parm, Data) {
    shape <- interval(parm[1], 0.01, 100)
    parm[1] <- shape

    bf <- rep(0.25, 4)
    Q  <- rep(1, 6)
    res <- Data$phylo_loglik(bf, Q, shape, Data$tree, Data$pdat)
    LL  <- res$logLik

    shape_prior <- dexp(shape, rate = 1, log = TRUE)
    LP <- LL + shape_prior

    yhat <- rep(LL / Data$N, Data$N)
    list(LP = LP, Dev = -2 * LL, Monitor = c(LP, LL),
         yhat = yhat, parm = parm)
}

IV_jc <- c(1)

# Fit JC model
fit_jc <- lucifer(
    Model_jc, Data_jc, IV_jc,
    Iterations = 10000, Status = 5000, Thinning = 5,
    Algorithm = "AMWG",
    Specs = list(B = NULL, n = 0, Periodicity = 50)
)
```

### HKY+\\\Gamma\\ model

``` r

Data_hky <- list(
    mon.names  = c("LP", "LL"),
    parm.names = c("alpha_A", "alpha_C", "alpha_G", "kappa", "shape"),
    tree = tree_rooted,
    pdat = pdat,
    N    = n_sites,
    phylo_loglik = phylo_loglik
)

Model_hky <- function(parm, Data) {
    # Base frequencies via softmax
    alpha <- c(parm[1:3], 0)
    alpha <- alpha - max(alpha)
    bf    <- exp(alpha) / sum(exp(alpha))

    # Transition/transversion ratio
    kappa <- interval(parm[4], 0.01, 200)
    parm[4] <- kappa

    # Gamma shape
    shape <- interval(parm[5], 0.01, 100)
    parm[5] <- shape

    # HKY: transitions (AG, CT) have rate kappa; transversions rate 1
    Q <- c(1, kappa, 1, 1, kappa, 1)

    res <- Data$phylo_loglik(bf, Q, shape, Data$tree, Data$pdat)
    LL  <- res$logLik

    # Priors
    pi_prior    <- sum(log(bf))  # Dirichlet(1,1,1,1) + Jacobian
    kappa_prior <- dexp(kappa, rate = 0.1, log = TRUE)
    shape_prior <- dexp(shape, rate = 1, log = TRUE)
    LP <- LL + pi_prior + kappa_prior + shape_prior

    yhat <- rep(LL / Data$N, Data$N)
    list(LP = LP, Dev = -2 * LL, Monitor = c(LP, LL),
         yhat = yhat, parm = parm)
}

IV_hky <- c(0, 0, 0, 2, 1)

# Fit HKY model
fit_hky <- lucifer(
    Model_hky, Data_hky, IV_hky,
    Iterations = 15000, Status = 5000, Thinning = 5,
    Algorithm = "AMWG",
    Specs = list(B = NULL, n = 0, Periodicity = 50)
)
```

### Model comparison via DIC

The deviance information criterion (DIC) provides a straightforward
comparison across models fitted to the same data. Lower DIC indicates
better predictive performance, with the effective number of parameters
\\p_D\\ penalizing model complexity. Note that SMC marginal likelihoods
would be the gold standard here, but
[`phangorn::pml()`](https://klausvigo.github.io/phangorn/reference/pml.html)
has a C-level memory management issue that causes segmentation faults
when called many thousands of times in rapid succession (as SMC requires
with its particle population). MCMC fits, which call
[`pml()`](https://klausvigo.github.io/phangorn/reference/pml.html)
sequentially with normal garbage-collection cycles, do not trigger this
problem.

``` r

# DIC comparison from MCMC fits
# DIC1 is a 3-element vector: (Dbar, pD, DIC)
dic_table <- data.frame(
    Model = c("JC69+G", "HKY+G", "GTR+G"),
    pD  = c(fit_jc$DIC1[2], fit_hky$DIC1[2], fit_gtr$DIC1[2]),
    DIC = c(fit_jc$DIC1[3], fit_hky$DIC1[3], fit_gtr$DIC1[3])
)
dic_table$Delta_DIC <- dic_table$DIC - min(dic_table$DIC)
dic_table <- dic_table[order(dic_table$DIC), ]

knitr::kable(dic_table, digits = 2, row.names = FALSE,
             caption = "DIC comparison of substitution models")
```

### Model comparison via LOO-PSIS

An alternative approach constructs per-site log-likelihood matrices from
each MCMC fit and applies LOO-PSIS. This is computationally cheaper than
SMC (no refitting) and provides per-site Pareto \\\hat{k}\\ diagnostics
that identify observations with outsized influence on the posterior.

``` r

# Helper: softmax transform for posterior alpha samples
softmax_bf <- function(alpha_vec) {
    alpha <- c(alpha_vec, 0)
    alpha <- alpha - max(alpha)
    exp(alpha) / sum(exp(alpha))
}

# Build per-site log-likelihood matrices for each model
build_jc_site_ll <- function(fit) {
    post <- as.matrix(fit$Posterior2)
    if (ncol(post) == 1) colnames(post) <- "shape"
    S <- nrow(post)
    ll_matrix <- matrix(NA_real_, nrow = n_sites, ncol = S)
    for (s in seq_len(S)) {
        res <- phylo_loglik(rep(0.25, 4), rep(1, 6), post[s, 1],
                            tree_rooted, pdat)
        ll_matrix[, s] <- res$site_ll
    }
    ll_matrix
}

build_hky_site_ll <- function(fit) {
    post <- fit$Posterior2; S <- nrow(post)
    ll_matrix <- matrix(NA_real_, nrow = n_sites, ncol = S)
    for (s in seq_len(S)) {
        bf    <- softmax_bf(post[s, 1:3])
        kappa <- post[s, 4]; shape <- post[s, 5]
        Q     <- c(1, kappa, 1, 1, kappa, 1)
        res   <- phylo_loglik(bf, Q, shape, tree_rooted, pdat)
        ll_matrix[, s] <- res$site_ll
    }
    ll_matrix
}

build_gtr_site_ll <- function(fit) {
    post <- fit$Posterior2; S <- nrow(post)
    ll_matrix <- matrix(NA_real_, nrow = n_sites, ncol = S)
    for (s in seq_len(S)) {
        bf    <- softmax_bf(post[s, 1:3])
        Q     <- c(post[s, 4:8], 1)
        shape <- post[s, 9]
        res   <- phylo_loglik(bf, Q, shape, tree_rooted, pdat)
        ll_matrix[, s] <- res$site_ll
    }
    ll_matrix
}

# Compute LOO-PSIS for each model
ll_mat_jc  <- build_jc_site_ll(fit_jc)
ll_mat_hky <- build_hky_site_ll(fit_hky)
ll_mat_gtr <- build_gtr_site_ll(fit_gtr)

loo_jc  <- LOO(ll_mat_jc)
loo_hky <- LOO(ll_mat_hky)
loo_gtr <- LOO(ll_mat_gtr)

# Compare all three
loo_comp <- loo_compare(list(JC = loo_jc, HKY = loo_hky, GTR = loo_gtr))
print(loo_comp)

# LOO diagnostics for best model
summary(loo_gtr)
plot(loo_gtr)
```

### Stacking weights

When no single model dominates, Bayesian stacking [\[17\]](#ref17) finds
the optimal combination of models that maximizes the leave-one-out
predictive density. This is more robust than selecting a single best
model, especially when models capture complementary aspects of the data.

``` r

weights <- stacking_weights(list(JC = loo_jc, HKY = loo_hky, GTR = loo_gtr))
print(weights)
```

## Example 3: strict versus relaxed molecular clock

This example demonstrates hierarchical Bayesian inference for molecular
clock models, one of the core problems in phylodynamics. We estimate
branch-specific evolutionary rates under both strict and relaxed clock
assumptions, using the HKY+\\\Gamma\\ substitution model for simplicity.

### Time-calibrated tree

We first construct an ultrametric (time-calibrated) tree using penalized
likelihood. The branch lengths of this tree represent time rather than
substitutions; the clock model then relates time to substitutions
through the rate parameter(s).

``` r

# Construct ultrametric tree from NJ tree
# Replace any zero-length branches (common in NJ trees) with a small value
# before calibrating, as chronopl requires strictly positive branch lengths
tree_for_cal <- tree_rooted
tree_for_cal$edge.length[tree_for_cal$edge.length <= 0] <- 1e-6

# Penalized likelihood calibration
tree_ultra <- ape::chronopl(tree_for_cal, lambda = 1, age.min = 1)

# Verify ultrametricity
cat("Is ultrametric:", ape::is.ultrametric(tree_ultra), "\n")
cat("N branches:", nrow(tree_ultra$edge), "\n")
n_branches <- nrow(tree_ultra$edge)

# Store time-branch lengths (will be scaled by rates)
time_lengths <- tree_ultra$edge.length
cat("Time branch lengths[1:5]:", round(time_lengths[1:5], 4), "\n")
```

### Strict clock model

Under the strict clock, all branches share a single evolutionary rate
\\r\\. The substitution branch length for branch \\i\\ is \\b_i = r
\cdot t_i\\.

``` r

Data_strict <- list(
    mon.names  = c("LP", "LL"),
    parm.names = c("rate", "kappa", "alpha_A", "alpha_C", "alpha_G", "shape"),
    tree = tree_ultra,
    pdat = pdat,
    time_lengths = time_lengths,
    N = n_sites,
    phylo_loglik = phylo_loglik
)

Model_strict <- function(parm, Data) {
    # Global clock rate
    rate <- interval(parm[1], 1e-6, 100)
    parm[1] <- rate

    # HKY parameters
    kappa <- interval(parm[2], 0.01, 200)
    parm[2] <- kappa

    # Base frequencies via softmax
    alpha <- c(parm[3:5], 0)
    alpha <- alpha - max(alpha)
    bf    <- exp(alpha) / sum(exp(alpha))

    shape <- interval(parm[6], 0.01, 100)
    parm[6] <- shape

    Q <- c(1, kappa, 1, 1, kappa, 1)

    # Scale time-tree by global rate
    scaled_tree <- Data$tree
    scaled_tree$edge.length <- Data$time_lengths * rate

    res <- Data$phylo_loglik(bf, Q, shape, scaled_tree, Data$pdat)
    LL  <- res$logLik

    # Priors
    rate_prior  <- dlnorm(rate, meanlog = 0, sdlog = 2, log = TRUE)
    kappa_prior <- dexp(kappa, rate = 0.1, log = TRUE)
    pi_prior    <- sum(log(bf))  # Dirichlet(1,1,1,1) + Jacobian
    shape_prior <- dexp(shape, rate = 1, log = TRUE)

    LP <- LL + rate_prior + kappa_prior + pi_prior + shape_prior

    yhat <- rep(LL / Data$N, Data$N)
    list(LP = LP, Dev = -2 * LL, Monitor = c(LP, LL),
         yhat = yhat, parm = parm)
}

IV_strict <- c(1, 5, 0, 0, 0, 0.5)

# Fit with AMWG (phangorn likelihood, no analytical gradients)
fit_strict <- lucifer(
    Model_strict, Data_strict, IV_strict,
    Iterations = 20000, Status = 5000, Thinning = 10,
    Algorithm = "AMWG",
    Specs = list(B = NULL, n = 0, Periodicity = 50)
)

print(fit_strict)
```

### Relaxed lognormal clock model

Under the relaxed clock, each branch has its own rate drawn from a
lognormal distribution with estimated hyperparameters. This is a
hierarchical model with \\n\_{\text{branches}} + 2\\ clock parameters
(the hyperparameters \\\mu_r, \sigma_r\\ plus the branch-specific rates)
in addition to the substitution model parameters.

We use a non-centered parameterization for the branch rates: instead of
sampling \\r_i\\ directly, we sample standardized scores \\z_i \sim
\mathcal{N}(0,1)\\ and recover rates as \\r_i = \exp(\mu_r + \sigma_r
z_i)\\. This reparameterization puts all branch parameters on a common
unit scale, which improves AMWG’s componentwise adaptation since the
step size for each \\z_i\\ operates on the same natural scale regardless
of the overall rate magnitude. The lognormal prior density and the
Jacobian of the exponential transformation cancel exactly, so the prior
on each \\z_i\\ reduces to a standard normal.

``` r

parm_names_relaxed <- c(
    "mu_r", "sigma_r",
    paste0("z_", seq_len(n_branches)),
    "kappa", "alpha_A", "alpha_C", "alpha_G", "shape"
)

Data_relaxed <- list(
    mon.names  = c("LP", "LL", "mean_rate"),
    parm.names = parm_names_relaxed,
    tree = tree_ultra,
    pdat = pdat,
    time_lengths = time_lengths,
    n_branches = n_branches,
    N = n_sites,
    phylo_loglik = phylo_loglik
)

Model_relaxed <- function(parm, Data) {
    nb <- Data$n_branches

    # Hyperparameters for branch rate distribution
    mu_r    <- parm[1]  # log-scale mean (unconstrained)
    sigma_r <- interval(parm[2], 0.01, 5)
    parm[2] <- sigma_r

    # Non-centered parameterization: z ~ N(0,1), rate = exp(mu_r + sigma_r * z)
    # The lognormal density and exponential Jacobian cancel, so the
    # prior on z is simply standard normal.
    z     <- parm[3:(nb + 2)]
    rates <- exp(mu_r + sigma_r * z)

    # HKY substitution parameters
    idx <- nb + 2
    kappa <- interval(parm[idx + 1], 0.01, 200)
    parm[idx + 1] <- kappa

    # Base frequencies via softmax
    alpha <- c(parm[(idx + 2):(idx + 4)], 0)
    alpha <- alpha - max(alpha)
    bf    <- exp(alpha) / sum(exp(alpha))

    shape <- interval(parm[idx + 5], 0.01, 100)
    parm[idx + 5] <- shape

    Q <- c(1, kappa, 1, 1, kappa, 1)

    # Scale time-tree by branch-specific rates
    scaled_tree <- Data$tree
    scaled_tree$edge.length <- Data$time_lengths * rates

    res <- Data$phylo_loglik(bf, Q, shape, scaled_tree, Data$pdat)
    LL  <- res$logLik

    # --- Priors ---
    z_prior     <- sum(dnorm(z, 0, 1, log = TRUE))
    mu_prior    <- dnorm(mu_r, mean = 0, sd = 2, log = TRUE)
    sigma_prior <- dhalfcauchy(sigma_r, scale = 1, log = TRUE)
    kappa_prior <- dexp(kappa, rate = 0.1, log = TRUE)
    pi_prior    <- sum(log(bf))  # Dirichlet(1,1,1,1) + Jacobian
    shape_prior <- dexp(shape, rate = 1, log = TRUE)

    LP <- LL + z_prior + mu_prior + sigma_prior +
          kappa_prior + pi_prior + shape_prior

    yhat <- rep(LL / Data$N, Data$N)
    mean_rate <- mean(rates)

    list(LP = LP, Dev = -2 * LL, Monitor = c(LP, LL, mean_rate),
         yhat = yhat, parm = parm)
}

# Initial values: mu_r from strict clock posterior, z = 0 (at the mean),
# substitution parameters from reasonable starting points.
strict_rate_mean <- mean(fit_strict$Posterior2[, 1])
IV_relaxed <- c(
    log(strict_rate_mean), 0.5,   # mu_r, sigma_r
    rep(0, n_branches),           # z-scores (all branches at mean rate)
    5, 0, 0, 0, 0.5              # kappa, alpha_A, alpha_C, alpha_G, shape
)

# AMWG evaluates the Model once per parameter per iteration (35 x 8K =
# 280K pml calls). The phylo_loglik wrapper calls gc() every 100 pml()
# evaluations, which prevents the C-level heap corruption that phangorn
# accumulates during intensive repeated evaluation. The non-centered
# parameterization avoids scale mismatch on the branch rates.
fit_relaxed <- lucifer(
    Model_relaxed, Data_relaxed, IV_relaxed,
    Iterations = 8000, Status = 4000, Thinning = 5,
    Algorithm = "AMWG",
    Specs = list(B = NULL, n = 0, Periodicity = 50),
    Chains = 1
)

print(fit_relaxed)
```

### Clock model comparison

``` r

# Compare strict vs relaxed via DIC.
# Note: SMC marginal likelihoods would be preferable but phangorn::pml()
# triggers C-level heap corruption under the intensive repeated
# evaluation that SMC requires with its particle population.
# The strict clock may show an anomalously large pD (effective number of
# parameters) when the single-rate assumption is a poor fit: the posterior
# mean rate gives a much better deviance than typical posterior draws,
# inflating pD = Dbar - D(theta_bar) far beyond the nominal parameter count.
dic_clock <- data.frame(
    Model = c("Strict clock", "Relaxed clock"),
    pD  = c(fit_strict$DIC1[2], fit_relaxed$DIC1[2]),
    DIC = c(fit_strict$DIC1[3], fit_relaxed$DIC1[3])
)
dic_clock$Delta_DIC <- dic_clock$DIC - min(dic_clock$DIC)

knitr::kable(dic_clock, digits = 2, row.names = FALSE,
             caption = "Clock model comparison via DIC")

# Interpret
delta <- abs(diff(dic_clock$DIC))
best <- dic_clock$Model[which.min(dic_clock$DIC)]
if (delta > 10) {
    cat(sprintf("Strong evidence for %s (Delta DIC = %.1f).\n", best, delta))
} else if (delta > 5) {
    cat(sprintf("Moderate evidence for %s (Delta DIC = %.1f).\n", best, delta))
} else {
    cat(sprintf("Models are comparable (Delta DIC = %.1f).\n", delta))
}
```

### Rate variation across branches

``` r

# Posterior distribution of rate variation (sigma_r)
post_relaxed <- fit_relaxed$Posterior2
mu_r_post    <- post_relaxed[, 1]
sigma_r_post <- post_relaxed[, 2]

cat("Rate heterogeneity (sigma_r):\n")
cat(sprintf("  Mean: %.4f\n", mean(sigma_r_post)))
cat(sprintf("  95%% CI: [%.4f, %.4f]\n",
            quantile(sigma_r_post, 0.025),
            quantile(sigma_r_post, 0.975)))

# Reconstruct branch rates from non-centered z-scores:
# rate_i = exp(mu_r + sigma_r * z_i)
z_cols <- 3:(n_branches + 2)
z_post <- post_relaxed[, z_cols]
rates_post <- exp(mu_r_post + sigma_r_post * z_post)

rate_means <- colMeans(rates_post)

rate_df <- data.frame(
    branch = seq_len(n_branches),
    mean_rate = rate_means,
    lower = apply(rates_post, 2, quantile, 0.025),
    upper = apply(rates_post, 2, quantile, 0.975)
)

ggplot2::ggplot(rate_df, ggplot2::aes(x = branch, y = mean_rate)) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper),
                           width = 0.3, alpha = 0.5) +
    ggplot2::geom_hline(yintercept = mean(rate_means),
                        linetype = "dashed", color = "steelblue") +
    ggplot2::labs(x = "Branch index", y = "Evolutionary rate",
                  title = "Branch-specific rate estimates (relaxed clock)") +
    ggplot2::theme_minimal()
```

## Example 4: coalescent demographic inference

This example estimates population demographic parameters from a
genealogy using the coalescent framework, a core application in
phylodynamics. Unlike the previous examples, no sequence likelihood is
needed; the coalescent provides a likelihood for the tree itself,
conditioned on the demographic model. We implement the coalescent
log-likelihood directly in R and fit three demographic models of
increasing flexibility.

### Coalescent log-likelihood

We define the coalescent log-likelihood for each demographic
parameterization. For a genealogy with \\n\\ tips, \\n - 1\\ coalescent
intervals with durations \\\Delta t_j\\ and lineage counts \\k_j\\, the
log-likelihood under a piecewise-constant demographic function is

\\\log L = \sum\_{j=1}^{n-1} \left\[ \log \binom{k_j}{2} - \log N_j -
\frac{\binom{k_j}{2} \Delta t_j}{N_j} \right\]\\

where \\N_j\\ is the effective population size during interval \\j\\.

``` r

#' Coalescent log-likelihood for piecewise-constant Ne
#' @param Ne Numeric vector of population sizes (one per interval or grouped)
#' @param intervals Numeric vector of interval durations (n-1 values)
#' @param lineages Integer vector of lineage counts (n-1 values)
#' @param groups Integer vector mapping intervals to Ne groups (optional)
#' @return Scalar log-likelihood
coal_loglik <- function(Ne, intervals, lineages, groups = NULL) {
    n_int <- length(intervals)
    if (is.null(groups)) groups <- seq_len(n_int)

    choose_k <- lineages * (lineages - 1) / 2
    ll <- 0
    for (j in seq_len(n_int)) {
        Ne_j <- Ne[groups[j]]
        ll <- ll + log(choose_k[j]) - log(Ne_j) -
              choose_k[j] * intervals[j] / Ne_j
    }
    ll
}
```

### Extracting coalescent intervals

We use the `bird.orders` tree from `ape`, an ultrametric phylogeny of 23
avian orders, to demonstrate coalescent inference. This tree is
ultrametric by construction, making it suitable for coalescent analysis.

``` r

# Load ultrametric tree
data(bird.orders, package = "ape")

# Extract coalescent intervals
ci <- ape::coalescent.intervals(bird.orders)
intervals <- ci$interval.length   # durations between coalescence events
lineages  <- ci$lineages          # number of lineages in each interval
n_tips    <- length(bird.orders$tip.label)
n_int     <- length(intervals)

cat("Tips:", n_tips, "\n")
cat("Coalescent intervals:", n_int, "\n")
cat("Total tree depth:", sum(intervals), "\n")
```

### Constant-size coalescent

The simplest demographic model assumes a constant effective population
size \\N_e\\ throughout the history of the sample. This model has one
free parameter.

``` r

Data_const <- list(
    mon.names  = c("LP"),
    parm.names = c("Ne"),
    intervals  = intervals,
    lineages   = lineages,
    n_int      = n_int,
    N          = n_int,  # "observations" = coalescent intervals
    coal_loglik = coal_loglik  # pass for callr visibility
)

Model_const <- function(parm, Data) {
    Ne <- interval(parm[1], 1e-6, 1e6)
    parm[1] <- Ne

    LL <- Data$coal_loglik(Ne, Data$intervals, Data$lineages,
                           groups = rep(1, Data$n_int))

    # Prior: log-uniform (Jeffreys prior on scale parameter)
    Ne_prior <- -log(Ne)

    LP <- LL + Ne_prior
    yhat <- rep(LL / Data$N, Data$N)

    list(LP = LP, Dev = -2 * LL, Monitor = LP,
         yhat = yhat, parm = parm)
}

# Analytical MLE (used for initial values and comparison)
choose_k <- lineages * (lineages - 1) / 2
Ne_mle <- sum(choose_k * intervals) / (n_tips - 1)

IV_const <- c(Ne_mle)

fit_const <- lucifer(
    Model_const, Data_const, IV_const,
    Iterations = 20000, Status = 5000, Thinning = 5,
    Algorithm = "NUTS",
    Specs = list(A = 5000, delta = 0.8, epsilon = NULL, Lmax = 10)
)

cat("MLE Ne:", Ne_mle, "\n")
cat("Posterior mean Ne:", mean(fit_const$Posterior2), "\n")
```

### Exponential growth coalescent

Under exponential growth, the population size at time \\t\\ (backward
from present) is \\N_e(t) = N_0 \exp(g t)\\, where \\N_0\\ is the
present-day size and \\g\\ is the growth rate. The coalescent
log-likelihood integral has a closed-form solution:
\\\int\_{t\_{j-1}}^{t_j} \binom{k}{2} / N_e(s) \\ ds = \binom{k}{2}
(e^{-g t\_{j-1}} - e^{-g t_j}) / (g N_0)\\ when \\g \neq 0\\. We
parameterize the sampler in terms of \\\log N_0\\ rather than \\N_0\\
directly; this reduces the strong correlation between population size
and growth rate that arises on the natural scale, improving NUTS
convergence.

``` r

#' Coalescent log-likelihood under exponential growth
coal_loglik_exp <- function(N0, g, intervals, lineages) {
    n_int <- length(intervals)
    choose_k <- lineages * (lineages - 1) / 2

    # Cumulative times (backward from present)
    cum_t <- c(0, cumsum(intervals))

    ll <- 0
    for (j in seq_len(n_int)) {
        t_start <- cum_t[j]
        t_end   <- cum_t[j + 1]
        Ne_t    <- N0 * exp(g * t_end)  # Ne at coalescence time

        ll <- ll + log(choose_k[j]) - log(Ne_t)

        if (abs(g) > 1e-10) {
            integral <- choose_k[j] *
                (exp(-g * t_start) - exp(-g * t_end)) / (g * N0)
        } else {
            integral <- choose_k[j] * intervals[j] / N0
        }
        ll <- ll - integral
    }
    ll
}

Data_exp <- list(
    mon.names  = c("LP"),
    parm.names = c("log_N0", "growth"),
    intervals  = intervals,
    lineages   = lineages,
    n_int      = n_int,
    N          = n_int,
    coal_loglik_exp = coal_loglik_exp  # pass for callr visibility
)

Model_exp <- function(parm, Data) {
    # Log-parameterization removes N0-g correlation and ensures positivity
    log_N0 <- parm[1]
    N0     <- exp(log_N0)
    g      <- parm[2]  # growth rate can be negative (decline)

    LL <- Data$coal_loglik_exp(N0, g, Data$intervals, Data$lineages)

    # Priors (both on unconstrained scale)
    log_N0_prior <- dnorm(log_N0, mean = 0, sd = 10, log = TRUE)
    g_prior      <- dnorm(g, mean = 0, sd = 2, log = TRUE)

    LP <- LL + log_N0_prior + g_prior
    yhat <- rep(LL / Data$N, Data$N)

    list(LP = LP, Dev = -2 * LL, Monitor = LP,
         yhat = yhat, parm = parm)
}

IV_exp <- c(log(Ne_mle), 0)

fit_exp <- lucifer(
    Model_exp, Data_exp, IV_exp,
    Iterations = 40000, Status = 10000, Thinning = 10,
    Algorithm = "NUTS",
    Specs = list(A = 10000, delta = 0.65, epsilon = NULL, Lmax = 20)
)

cat("Posterior N0:", mean(exp(fit_exp$Posterior2[, 1])), "\n")
cat("Posterior growth rate:", mean(fit_exp$Posterior2[, 2]), "\n")
```

### Bayesian skyline

The Bayesian skyline model divides the coalescent history into \\M\\
groups of consecutive intervals, each with its own population size
parameter. This non-parametric approach can recover complex demographic
trajectories. We use \\M = 5\\ groups with a random walk prior on the
log-population sizes to impose smoothness.

``` r

M <- 5  # number of population size groups

# Assign intervals to groups (approximately equal number per group)
ints_per_group <- ceiling(n_int / M)
groups <- rep(seq_len(M), each = ints_per_group)[seq_len(n_int)]

# Group midpoint times for plotting
cum_t <- c(0, cumsum(intervals))
group_times <- tapply(
    (cum_t[seq_len(n_int)] + cum_t[seq_len(n_int) + 1]) / 2,
    groups, mean
)

Data_sky <- list(
    mon.names  = c("LP"),
    parm.names = paste0("log_Ne_", seq_len(M)),
    intervals  = intervals,
    lineages   = lineages,
    groups     = groups,
    n_int      = n_int,
    M          = M,
    group_times = as.numeric(group_times),
    N          = n_int,
    coal_loglik = coal_loglik  # pass for callr visibility
)

Model_sky <- function(parm, Data) {
    # Log-population sizes (unconstrained)
    log_Ne <- parm[seq_len(Data$M)]
    Ne <- exp(log_Ne)

    LL <- Data$coal_loglik(Ne, Data$intervals, Data$lineages, Data$groups)

    # Random walk prior on log-Ne (smoothness)
    rw_prior <- 0
    for (m in 2:Data$M) {
        rw_prior <- rw_prior +
            dnorm(log_Ne[m], mean = log_Ne[m - 1], sd = 1, log = TRUE)
    }
    # Broad prior on first log-Ne
    first_prior <- dnorm(log_Ne[1], mean = 0, sd = 5, log = TRUE)

    LP <- LL + rw_prior + first_prior
    yhat <- rep(LL / Data$N, Data$N)

    list(LP = LP, Dev = -2 * LL, Monitor = LP,
         yhat = yhat, parm = parm)
}

IV_sky <- rep(log(Ne_mle), M)  # start at log(MLE) for each group

fit_sky <- lucifer(
    Model_sky, Data_sky, IV_sky,
    Iterations = 30000, Status = 5000, Thinning = 10,
    Algorithm = "NUTS",
    Specs = list(A = 5000, delta = 0.8, epsilon = NULL, Lmax = 10)
)
```

### Skyline plot

The Bayesian skyline plot reconstructs effective population size through
time with uncertainty quantification, directly comparable to the output
from BEAST’s Bayesian skyline analysis.

``` r

post_sky <- fit_sky$Posterior2

# Compute Ne posteriors per group
sky_df <- data.frame(
    time = Data_sky$group_times,
    Ne_median = apply(exp(post_sky), 2, median),
    Ne_lower  = apply(exp(post_sky), 2, quantile, 0.025),
    Ne_upper  = apply(exp(post_sky), 2, quantile, 0.975)
)

ggplot2::ggplot(sky_df, ggplot2::aes(x = time)) +
    ggplot2::geom_ribbon(
        ggplot2::aes(ymin = Ne_lower, ymax = Ne_upper),
        fill = "steelblue", alpha = 0.3
    ) +
    ggplot2::geom_line(ggplot2::aes(y = Ne_median),
                       color = "steelblue", linewidth = 1.2) +
    ggplot2::geom_point(ggplot2::aes(y = Ne_median),
                        color = "steelblue", size = 3) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_x_reverse() +
    ggplot2::labs(
        x = "Time (coalescent units, past to present)",
        y = expression(paste("Effective population size  ", N[e])),
        title = "Bayesian skyline reconstruction"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 14),
        axis.title = ggplot2::element_text(size = 12)
    )
```

### Demographic model comparison

``` r

# Compare demographic models via SMC
smc_const <- SMC(Model_const, Data_const, IV_const,
                 N.particles = 500, Rejuvenation.steps = 5)
smc_exp   <- SMC(Model_exp, Data_exp, IV_exp,
                 N.particles = 500, Rejuvenation.steps = 5)
smc_sky   <- SMC(Model_sky, Data_sky, IV_sky,
                 N.particles = 500, Rejuvenation.steps = 5)

cat("Log marginal likelihoods:\n")
cat(sprintf("  Constant:    %.2f\n", smc_const$Log.Marginal.Likelihood))
cat(sprintf("  Exponential: %.2f\n", smc_exp$Log.Marginal.Likelihood))
cat(sprintf("  Skyline:     %.2f\n", smc_sky$Log.Marginal.Likelihood))
```

## Example 5: multimodality detection and diagnostics

Gao et al. [\[3\]](#ref3) demonstrated that phylodynamic posterior
landscapes are inherently and extensively multimodal: their analysis of
15 large viral datasets showed that MCMC chains routinely become trapped
in separate peaks of tree space, with deep valleys that require
exceedingly rare SPR moves to cross. Standard continuous-parameter
diagnostics (ESS on substitution rates, clock parameters) often fail to
detect these problems, while tree-specific diagnostics (parsimony-score
ESS, clade-specific \\\hat{R}\\) reveal them.

This finding has a direct implication for the conditional-inference
approach used in this vignette: if we condition on the wrong tree
topology, our parameter estimates may be biased. The solution is
transparency. Rather than relying on a single tree from a joint sampler
whose topology mixing may be inadequate, we can explicitly condition on
multiple candidate trees, fit the same model to each, and assess whether
the parameter estimates are sensitive to the tree choice. If they are,
we know that tree uncertainty matters for those parameters; if they are
not, we have evidence of robustness.

### Setup: two candidate trees

We construct two different trees from the woodmouse data using different
methods, then fit the GTR+\\\Gamma\\ model (from Example 1) to each.

``` r

# Tree 1: Neighbor-joining (already constructed)
tree1 <- tree_rooted

# Tree 2: Maximum parsimony (topology), then add branch lengths via pml
tree2_unrooted <- phangorn::pratchet(pdat, trace = 0, all = FALSE)
tree2_unrooted <- ape::compute.brlen(tree2_unrooted)  # add initial edge weights
tree2_fit <- phangorn::pml(tree2_unrooted, pdat)      # optimise branch lengths
tree2_unrooted <- tree2_fit$tree
tree2 <- ape::root(tree2_unrooted, outgroup = "No305", resolve.root = TRUE)

# Compare topologies
rf_dist <- phangorn::RF.dist(tree1, tree2)
cat("Robinson-Foulds distance between trees:", rf_dist, "\n")
```

### Parallel chains with dispersed initialization

For each tree, we run 4 parallel chains from dispersed starting points.
This is the standard approach for detecting multimodality: if chains
converge to different regions of parameter space, the posterior is
likely multimodal.

``` r

# Generate dispersed initial values as a matrix (one row per chain)
set.seed(42)
n_chains <- 4
init_matrix <- t(replicate(n_chains, {
    c(
        rnorm(3, 0, 1),               # softmax alphas (dispersed)
        rexp(5, rate = 0.5),          # GTR rates
        rexp(1, rate = 0.5)           # gamma shape
    )
}))

# Fit on tree 1 with 4 parallel chains
Data_tree1 <- Data  # already uses tree_rooted = tree1

fit_tree1 <- lucifer(
    Model, Data_tree1, init_matrix,
    Iterations = 20000, Thinning = 10,
    Algorithm = "AMWG",
    Specs = list(B = NULL, n = 0, Periodicity = 50),
    Chains = n_chains, CPUs = n_chains
)

# Fit on tree 2
Data_tree2 <- Data
Data_tree2$tree <- tree2

fit_tree2 <- lucifer(
    Model, Data_tree2, init_matrix,
    Iterations = 20000, Thinning = 10,
    Algorithm = "AMWG",
    Specs = list(B = NULL, n = 0, Periodicity = 50),
    Chains = n_chains, CPUs = n_chains
)
```

### Convergence diagnostics across chains

``` r

# Extract per-chain posterior matrices
posts_t1 <- lapply(fit_tree1$Chains, function(ch) ch$Posterior2)
posts_t2 <- lapply(fit_tree2$Chains, function(ch) ch$Posterior2)

# Compute R-hat across chains for each parameter
compute_multichain_rhat <- function(posts, parm_names) {
    # Truncate chains to equal length (thinning can produce slight mismatches)
    min_n <- min(sapply(posts, nrow))
    posts <- lapply(posts, function(x) x[seq_len(min_n), , drop = FALSE])

    n_parms <- ncol(posts[[1]])
    rhat_vals <- numeric(n_parms)
    ess_bulk_vals <- numeric(n_parms)

    for (p in seq_len(n_parms)) {
        chains <- do.call(cbind, lapply(posts, function(x) x[, p]))
        rhat_vals[p] <- Rhat(chains)
        ess_bulk_vals[p] <- ESS.bulk(as.vector(chains))
    }

    data.frame(
        Parameter = parm_names,
        Rhat = round(rhat_vals, 4),
        ESS.bulk = round(ess_bulk_vals, 0)
    )
}

cat("--- Diagnostics on Tree 1 (NJ) ---\n")
diag_t1 <- compute_multichain_rhat(posts_t1, Data$parm.names)
print(diag_t1)

cat("\n--- Diagnostics on Tree 2 (MP) ---\n")
diag_t2 <- compute_multichain_rhat(posts_t2, Data$parm.names)
print(diag_t2)
```

### Sensitivity to tree topology

We compare the posterior distributions of key parameters across the two
trees. If the posteriors overlap substantially, the parameter estimates
are robust to tree choice; if they diverge, tree uncertainty dominates.

``` r

# Combine posteriors from all chains for each tree
all_post_t1 <- do.call(rbind, posts_t1)
all_post_t2 <- do.call(rbind, posts_t2)

# Helper: transform softmax alphas to pi_A for a posterior matrix
softmax_pi_A <- function(post_mat) {
    apply(post_mat[, 1:3, drop = FALSE], 1, function(a) {
        alpha <- c(a, 0); alpha <- alpha - max(alpha)
        exp(alpha[1]) / sum(exp(alpha))
    })
}

# Compare posterior distributions for each parameter
# For alphas, report the transformed frequencies instead
bf_t1 <- t(apply(all_post_t1[, 1:3], 1, function(a) {
    alpha <- c(a, 0); alpha <- alpha - max(alpha)
    exp(alpha) / sum(exp(alpha))
}))
bf_t2 <- t(apply(all_post_t2[, 1:3], 1, function(a) {
    alpha <- c(a, 0); alpha <- alpha - max(alpha)
    exp(alpha) / sum(exp(alpha))
}))

# Build comparison table with transformed frequencies
parm_labels <- c("pi_A", "pi_C", "pi_G",
                 Data$parm.names[4:9])
combined_t1 <- cbind(bf_t1[, 1:3], all_post_t1[, 4:9])
combined_t2 <- cbind(bf_t2[, 1:3], all_post_t2[, 4:9])

par_compare <- data.frame(
    Parameter = parm_labels,
    Tree1_mean = colMeans(combined_t1),
    Tree1_sd   = apply(combined_t1, 2, sd),
    Tree2_mean = colMeans(combined_t2),
    Tree2_sd   = apply(combined_t2, 2, sd)
)
par_compare$Overlap <- with(par_compare, {
    sigma_pooled <- sqrt((Tree1_sd^2 + Tree2_sd^2) / 2)
    diff_means <- abs(Tree1_mean - Tree2_mean)
    exp(-0.25 * diff_means^2 / sigma_pooled^2)
})

knitr::kable(par_compare, digits = 4,
             caption = "Parameter estimates conditional on two alternative trees")

# Density comparison plot for key parameters
plot_data <- rbind(
    data.frame(tree = "NJ", param = "q_AG",
               value = all_post_t1[, 5]),
    data.frame(tree = "MP", param = "q_AG",
               value = all_post_t2[, 5]),
    data.frame(tree = "NJ", param = "shape",
               value = all_post_t1[, 9]),
    data.frame(tree = "MP", param = "shape",
               value = all_post_t2[, 9]),
    data.frame(tree = "NJ", param = "pi_A",
               value = softmax_pi_A(all_post_t1)),
    data.frame(tree = "MP", param = "pi_A",
               value = softmax_pi_A(all_post_t2))
)

ggplot2::ggplot(plot_data, ggplot2::aes(x = value, fill = tree)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::facet_wrap(~ param, scales = "free", ncol = 1) +
    ggplot2::scale_fill_manual(values = c("NJ" = "steelblue",
                                          "MP" = "coral")) +
    ggplot2::labs(x = "Parameter value", y = "Density",
                  fill = "Tree",
                  title = "Posterior sensitivity to tree topology") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 14),
        strip.text = ggplot2::element_text(face = "bold", size = 12)
    )
```

### Prior sensitivity via power-scaling

lucifer’s
[`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md)
function provides additional diagnostics for detecting prior-likelihood
conflict through power-scaling sensitivity analysis [\[21\]](#ref21).
The idea is simple: if inflating the prior (\\\alpha \> 1\\) or
deflating it (\\\alpha \< 1\\) substantially changes the posterior, the
prior is influential and may dominate the likelihood. When the prior
strongly pulls estimates away from the likelihood, it may indicate model
misspecification or prior choice problems. In the phylogenetic context,
prior sensitivity of substitution rate parameters can reveal whether the
data contain sufficient information to estimate rates independently of
prior assumptions.

``` r

# RobustBayes power-scaling on the GTR model (tree 1):
# inflate/deflate prior (alpha < 1 weakens, alpha > 1 strengthens)
# to assess how sensitive posteriors are to prior choice.
# Uses PSIS importance reweighting, no refitting required.
rb <- RobustBayes(
    fit_tree1, Model = Model, Data = Data,
    modules = "power",
    alpha = seq(0.5, 1.5, by = 0.1)
)
summary(rb)
plot(rb)
```

### Connecting to Gao et al.

The diagnostics demonstrated here, parallel chains with \\\hat{R}\\,
density comparisons across conditioning trees, and prior sensitivity
analysis, address the core concern raised by Gao et al. [\[3\]](#ref3).
Their key findings were that (1) tree-topology multimodality is
widespread and biologically driven, (2) standard continuous-parameter
diagnostics miss it, and (3) the impact on downstream estimates (rates,
divergence times, demographic trajectories) can be substantial.

lucifer’s approach of explicit conditioning offers three advantages.
First, by fitting the same model to different candidate trees and
comparing parameter estimates, we directly quantify the sensitivity of
downstream inference to tree choice, something BEAST’s joint sampler
obscures when chains fail to mix across tree peaks. Second, lucifer’s
modern rank-normalized \\\hat{R}\\ and bulk/tail ESS diagnostics are
more sensitive to mixing problems than the classical ESS metrics used in
Tracer. Third, the availability of multiple MCMC algorithms allows the
analyst to test whether apparent convergence is algorithm-specific: if
NUTS and AIES give different results for the same model, this flags
potential problems that a single-algorithm approach would miss.

## Discussion

This vignette demonstrated that lucifer, a general-purpose Bayesian
inference engine, can address core phylodynamic problems by conditioning
on fixed tree topologies and leveraging its algorithm diversity, modern
diagnostics, and model comparison tools.

**What lucifer provides.** For continuous-parameter phylogenetic
inference, lucifer offers over 40 MCMC algorithms, ranging from
gradient-based methods (NUTS, HMC, MALA) through adaptive random-walk
methods (AMWG, DRAM, AM) to ensemble methods (AIES, DEMC, t-walk), each
with different strengths for different posterior geometries. The
Crucible pipeline automates algorithm selection, fitting multiple
methods and comparing their performance. Sequential Monte Carlo provides
unbiased marginal likelihood estimates for Bayes factor computation,
paralleling BEAST’s stepping-stone and path-sampling approaches but with
the additional capability of adaptive tempering. LOO-PSIS and stacking
enable model comparison and combination with per-observation diagnostics
that BEAST does not offer. Modern convergence diagnostics
(rank-normalized \\\hat{R}\\, bulk and tail ESS) are more sensitive to
mixing problems than the classical metrics available in Tracer
[\[18\]](#ref18).

**What lucifer does not provide.** lucifer has no tree topology sampler.
It cannot explore tree space, and it cannot jointly estimate trees and
parameters. This means it cannot replace BEAST or RevBayes
[\[19\]](#ref19) for problems where tree uncertainty is the primary
quantity of interest. It also lacks specialized phylogenetic features
such as tip-dating calibration, structured coalescent models,
birth-death epidemiological models, and BEAGLE GPU-accelerated
likelihood computation.

**The case for conditional inference.** As Gao et al. [\[3\]](#ref3)
showed, the joint sampling approach is not without its own problems.
When tree-topology mixing fails, as it does for 12 of 15 large
phylodynamic datasets in their analysis, the joint posterior is
effectively conditioned on whatever tree peak the MCMC chain happened to
find, without the analyst knowing it. Explicit conditioning on candidate
trees, with transparent sensitivity analysis, trades the aspiration of
full posterior integration for the honesty of knowing exactly what you
are conditioning on. This is particularly valuable in applied settings
where downstream conclusions (e.g., the geographic origin of an
epidemic, the timing of a population bottleneck) depend on whether
tree-topology uncertainty has been adequately captured.

**Future directions.** Three extensions could expand lucifer’s utility
in this domain. First, the product space sampler already implemented in
lucifer could be adapted for discrete model selection over a finite set
of candidate trees, enabling a form of topology averaging without
requiring SPR/NNI proposals. Second, lucifer could serve as a
post-processing engine for BEAST output, taking a posterior sample of
trees and estimating downstream parameters conditional on each tree,
then averaging across the tree posterior. Third, the SMC framework could
be extended with tree-aware mutation kernels, although this would
require substantial new infrastructure.

## References

**\[1\]** Suchard, M. A., Lemey, P., Baele, G., Ayres, D. L., Drummond,
A. J. & Rambaut, A. (2018). Bayesian phylogenetic and phylodynamic data
integration using BEAST 1.10. *Virus Evolution*, 4(1), vey016.
[doi:10.1093/ve/vey016](https://doi.org/10.1093/ve/vey016)

**\[2\]** Baele, G., Gill, M. S., Lemey, P. & Suchard, M. A. (2025).
BEAST X for Bayesian phylogenetic, phylogeographic and phylodynamic
inference. *Nature Methods*, 22, 1-4.
[doi:10.1038/s41592-024-02563-1](https://doi.org/10.1038/s41592-024-02563-1)

**\[3\]** Gao, J., Brusselmans, M., Carvalho, L. M., Suchard, M. A.,
Baele, G. & Matsen IV, F. A. (2026). Biological causes and impacts of
rugged tree landscapes in phylodynamic inference. *PNAS*, 123(2),
e2510938123. [DOI:
10.1073/pnas.2510938123](https://doi.org/10.1073/pnas.2510938123).

**\[4\]** Paradis, E. & Schliep, K. (2019). ape 5.0: an environment for
modern phylogenetics and evolutionary analyses in R. *Bioinformatics*,
35(3), 526-528.
[doi:10.1093/bioinformatics/bty633](https://doi.org/10.1093/bioinformatics/bty633)

**\[5\]** Schliep, K. P. (2011). phangorn: phylogenetic analysis in R.
*Bioinformatics*, 27(4), 592-593.
[doi:10.1093/bioinformatics/btq706](https://doi.org/10.1093/bioinformatics/btq706)

**\[6\]** Tavare, S. (1986). Some probabilistic and statistical problems
in the analysis of DNA sequences. *Lectures on Mathematics in the Life
Sciences*, 17, 57-86.

**\[7\]** Jukes, T. H. & Cantor, C. R. (1969). Evolution of protein
molecules. In H. N. Munro (Ed.), *Mammalian Protein Metabolism* (Vol.
III, pp. 21-132). Academic Press.

**\[8\]** Hasegawa, M., Kishino, H. & Yano, T. (1985). Dating of the
human-ape splitting by a molecular clock of mitochondrial DNA. *Journal
of Molecular Evolution*, 22(2), 160-174.
[doi:10.1007/BF02101694](https://doi.org/10.1007/BF02101694)

**\[9\]** Yang, Z. (1994). Maximum likelihood phylogenetic estimation
from DNA sequences with variable rates over sites: approximate methods.
*Journal of Molecular Evolution*, 39(3), 306-314.
[doi:10.1007/BF00160154](https://doi.org/10.1007/BF00160154)

**\[10\]** Felsenstein, J. (1981). Evolutionary trees from DNA
sequences: a maximum likelihood approach. *Journal of Molecular
Evolution*, 17(6), 368-376.
[doi:10.1007/BF01734359](https://doi.org/10.1007/BF01734359)

**\[11\]** Zuckerkandl, E. & Pauling, L. (1965). Evolutionary divergence
and convergence in proteins. In V. Bryson & H. J. Vogel (Eds.),
*Evolving Genes and Proteins* (pp. 97-166). Academic Press.

**\[12\]** Drummond, A. J., Ho, S. Y. W., Phillips, M. J. & Rambaut, A.
(2006). Relaxed phylogenetics and dating with confidence. *PLoS
Biology*, 4(5), e88.
[doi:10.1371/journal.pbio.0040088](https://doi.org/10.1371/journal.pbio.0040088)

**\[13\]** Kingman, J. F. C. (1982). The coalescent. *Stochastic
Processes and their Applications*, 13(3), 235-248.
[doi:10.1016/0304-4149(82)90011-4](https://doi.org/10.1016/0304-4149(82)90011-4)

**\[14\]** Drummond, A. J., Rambaut, A., Shapiro, B. & Pybus, O. G.
(2005). Bayesian coalescent inference of past population dynamics from
molecular sequences. *Molecular Biology and Evolution*, 22(5),
1185-1192.
[doi:10.1093/molbev/msi103](https://doi.org/10.1093/molbev/msi103)

**\[15\]** Michaux, J. R., Magnanou, E., Paradis, E., Nieberding, C. &
Libois, R. (2003). Mitochondrial phylogeography of the woodmouse
(*Apodemus sylvaticus*) in the Western Palearctic region. *Molecular
Ecology*, 12(3), 685-697.
[doi:10.1046/j.1365-294X.2003.01752.x](https://doi.org/10.1046/j.1365-294X.2003.01752.x)

**\[16\]** Hoffman, M. D. & Gelman, A. (2014). The No-U-Turn Sampler:
adaptively setting path lengths in Hamiltonian Monte Carlo. *Journal of
Machine Learning Research*, 15(1), 1593-1623.

**\[17\]** Yao, Y., Vehtari, A., Simpson, D. & Gelman, A. (2018). Using
stacking to average Bayesian predictive distributions. *Bayesian
Analysis*, 13(3), 917-1007.
[doi:10.1214/17-BA1091](https://doi.org/10.1214/17-BA1091)

**\[18\]** Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. &
Burkner, P.-C. (2021). Rank-normalization, folding, and localization: an
improved \\\hat{R}\\ for assessing convergence of MCMC (with
discussion). *Bayesian Analysis*, 16(2), 667-718.
[doi:10.1214/20-BA1221](https://doi.org/10.1214/20-BA1221)

**\[19\]** Hoehna, S., Landis, M. J., Heath, T. A., Boussau, B.,
Lartillot, N., Moore, B. R., Huelsenbeck, J. P. & Ronquist, F. (2016).
RevBayes: Bayesian phylogenetic inference using graphical models and an
interactive model-specification language. *Systematic Biology*, 65(4),
726-736.
[doi:10.1093/sysbio/syw021](https://doi.org/10.1093/sysbio/syw021)

**\[20\]** Roberts, G. O. & Rosenthal, J. S. (2009). Examples of
adaptive MCMC. *Journal of Computational and Graphical Statistics*,
18(2), 349-367.
[doi:10.1198/jcgs.2009.06134](https://doi.org/10.1198/jcgs.2009.06134)

**\[21\]** Kallioinen, N., Paananen, T., Burkner, P.-C. & Vehtari, A.
(2024). Detecting and diagnosing prior and likelihood sensitivity with
power-scaling. *Statistics and Computing*, 34, 57.
[doi:10.1007/s11222-023-10366-5](https://doi.org/10.1007/s11222-023-10366-5)
