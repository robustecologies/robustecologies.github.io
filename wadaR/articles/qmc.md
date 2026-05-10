# Quadratic Monte Carlo ensemble methods in lucifer

``` r

library(lucifer)
library(ggplot2)
set.seed(666)
```

  

## Introduction

Ensemble Markov chain Monte Carlo methods maintain a population of
interacting walkers that collectively explore a target distribution by
constructing proposals from the relative positions of other walkers. The
affine-invariant ensemble sampler (AIES) introduced by Goodman and Weare
[\[1\]](#ref1) proposes new positions by stretching the target walker
along the line connecting it to a randomly chosen companion, and this
linear move has become the default ensemble proposal through the widely
adopted **emcee** implementation [\[2\]](#ref2). Despite its simplicity
and affine invariance, the AIES stretch move interpolates linearly
between just two points; when the posterior concentrates along a curved,
narrow manifold, such as the banana-shaped distributions common in
hierarchical models or nonlinear regression, the linear proposal
overshoots the high-density region and acceptance rates drop sharply.

Militzer [\[3\]](#ref3) addressed this limitation by replacing the
linear interpolation with Lagrange polynomial interpolation through
three or more walkers, constructing proposals that follow the curvature
of the fitness landscape. A second-order (quadratic) interpolant through
three walkers traces a parabolic arc through the ensemble, and the
resulting Quadratic Monte Carlo (QMC) move can match banana-shaped
posteriors far better than a straight line. A subsequent paper
[\[4\]](#ref4) extended this idea with five additional move types that
combine higher-order polynomials, energy-directed proposals,
simplex-averaged guides, and symmetric random walks, yielding a family
of eight ensemble methods that cover a wide range of posterior
geometries.

The **lucifer** package implements all eight of Militzer’s ensemble
moves as MCMC algorithms backed by a C++ engine. Each method requires a
`Specs` list specifying the ensemble size `Nc`, a scale parameter `a`,
and method-specific options; the C++ backend handles walker
initialization, Lagrange weight computation, Jacobian correction, and
Metropolis-Hastings acceptance. Only walker 0 is recorded in the output
chain, matching the standard `demonoid` output interface, so all
existing diagnostic and plotting infrastructure applies directly.

This vignette develops the mathematical foundations of the QMC family in
detail, defines four benchmark distributions with known analytical
properties, runs all eight methods (plus AIES as baseline) on each
benchmark, and compares their performance through acceptance rates,
posterior recovery, effective sample size, and diagnostic plots. The
final section distills these results into practical guidance for
choosing among the methods.

  

## Mathematical foundations

### Lagrange polynomial interpolation

The central idea underlying the QMC family is to construct proposals by
polynomial interpolation through a set of walker positions. Given
\\n+1\\ points \\(t_0, \mathbf{r}\_0), (t_1, \mathbf{r}\_1), \ldots,
(t_n, \mathbf{r}\_n)\\ where each \\t_m\\ is a scalar knot value and
each \\\mathbf{r}\_m \in \mathbb{R}^N\\ is a walker position in
parameter space, the unique polynomial of degree \\n\\ passing through
all \\n+1\\ points is the Lagrange interpolant

\\ \mathbf{R}(t) = \sum\_{m=0}^{n} L_m(t) \\ \mathbf{r}\_m, \qquad
L_m(t) = \prod\_{\substack{j=0 \\ j \neq m}}^{n} \frac{t - t_j}{t_m -
t_j}. \\

Each basis polynomial \\L_m(t)\\ equals unity at \\t = t_m\\ and
vanishes at every other knot, so \\\mathbf{R}(t_m) = \mathbf{r}\_m\\
exactly. The Lagrange weights \\w_m = L_m(t')\\ evaluated at any
particular \\t'\\ satisfy the partition-of-unity property
\\\sum\_{m=0}^{n} w_m = 1\\, which follows from the fact that the
interpolant of the constant function \\f \equiv 1\\ must equal 1
everywhere. This means the proposed position \\\mathbf{R}(t') = \sum_m
w_m \mathbf{r}\_m\\ is an affine combination of the walker positions
(not necessarily convex, since individual weights can be negative),
ensuring that the proposal lies in the affine hull of the walkers and
inherits the affine invariance of the ensemble.

For the second-order case (QMC), three walkers define a quadratic arc;
the target walker \\i\\ receives a random knot \\t_i\\, while two guide
walkers receive fixed knots \\t_j = -1\\ and \\t_k = +1\\. The explicit
weights are

\\ w_i = \frac{(t' - t_j)(t' - t_k)}{(t_i - t_j)(t_i - t_k)}, \quad w_j
= \frac{(t' - t_i)(t' - t_k)}{(t_j - t_i)(t_j - t_k)}, \quad w_k =
\frac{(t' - t_i)(t' - t_j)}{(t_k - t_i)(t_k - t_j)}. \\

Third-order interpolation (QMC3) uses four walkers with knots at \\t_j =
-1\\, \\t_k = 0\\, \\t_l = +1\\, producing a cubic arc through four
points. The general \\n\\th-order case (QMCN) places \\n\\ guide walkers
at equispaced knots on \\\[-1, +1\]\\ and evaluates the Lagrange basis
of degree \\n\\ at a proposed knot \\t'\\. Higher-order interpolation
can capture more complex curvature but requires proportionally more
guide walkers and introduces more extreme Lagrange weights for knots far
from the interpolation nodes.

### Jacobian correction and detailed balance

The QMC proposal maps the old position \\\mathbf{r}\_i\\ to a new
position \\\mathbf{R}(t') = w_i \mathbf{r}\_i + \sum\_{m \neq i} w_m
\mathbf{r}\_m\\. Since the guide positions \\\mathbf{r}\_m\\ for \\m
\neq i\\ are fixed during this step, the transformation is a linear map
from \\\mathbf{r}\_i\\ to \\\mathbf{R}\\: each coordinate of
\\\mathbf{R}\\ depends on the corresponding coordinate of
\\\mathbf{r}\_i\\ through a multiplicative factor \\w_i\\ plus a
constant offset from the guide contributions. The Jacobian of this
affine map is therefore \\\|w_i\|^N\\, where \\N\\ is the parameter
dimension, and the acceptance probability that satisfies detailed
balance is

\\ \alpha = \min\left\[1, \\
\frac{\pi(\mathbf{R}(t'))}{\pi(\mathbf{r}\_i)} \\ \|w_i\|^N \right\]. \\

The Jacobian correction \\\|w_i\|^N\\ plays an analogous role to the
\\z^{N-1}\\ factor in AIES, penalizing proposals that deform the
parameter space too aggressively. When \\w_i \approx 1\\ (meaning the
proposed point is close to the original walker position with minor
contributions from the guides), the penalty vanishes; when \\w_i\\ is
large or small (meaning the proposal is dominated by guide walkers), the
Jacobian penalizes acceptance. Because the penalty scales exponentially
with dimension \\N\\, proposals with \\\|w_i\|\\ far from unity are
increasingly rejected in higher dimensions. This is why the scale
parameter `a`, which controls the spread of the knot distribution and
hence the distribution of \\w_i\\, must be reduced as \\N\\ grows.

### Knot sampling distributions

The knot assigned to the target walker is drawn from a distribution
\\P_S(t)\\ on \\\mathbb{R}\\, and the same distribution is used to
sample the proposal knot \\t'\\. The **lucifer** implementation offers
two choices controlled by the `Gaussian` flag in `Specs`:

The **uniform** distribution \\t \sim \text{Unif}\[-a, a\]\\ spreads
probability mass evenly across the knot range. When \\a\\ is moderate,
the Lagrange weight \\w_i\\ stays near unity most of the time, yielding
mild Jacobian penalties and reasonable acceptance rates. Uniform
sampling is recommended for low to moderate dimensions (\\N \leq 10\\)
and provides a simple, interpretable control over proposal
aggressiveness.

The **Gaussian** distribution \\t \sim \mathcal{N}(0, a^2)\\
concentrates mass near \\t = 0\\ and allows occasional large excursions.
For high-dimensional problems (\\N \> 10\\), Gaussian sampling naturally
favors identity-like proposals (where \\w_i \approx 1\\) while still
permitting exploratory moves; Militzer [\[4\]](#ref4) reports that this
produces more stable acceptance rates across dimensions compared to
uniform sampling.

### The stretch distribution for affine methods

The affine methods MAMC and SAMC use the stretch distribution from
Goodman and Weare [\[1\]](#ref1) rather than Lagrange interpolation. The
stretch factor \\z\\ is drawn from the density \\g(z) \propto
1/\sqrt{z}\\ supported on \\\[1/a, a\]\\ where \\a \> 1\\. This
distribution has the property that \\z\\ and \\1/z\\ have the same
density (up to the Jacobian of inversion), which is essential for
detailed balance. The mean of \\g\\ is not 1, so stretch moves produce a
mixture of compression (\\z \< 1\\) and expansion (\\z \> 1\\) relative
to the reference point. The standard recommendation \\a = 2\\ gives a
stretch range of \\\[0.5, 2\]\\ with most mass near \\z = 1\\, producing
a good balance between exploration and acceptance. The acceptance
probability for affine proposals with stretch factor \\z\\ in \\N\\
dimensions is

\\ \alpha = \min\left\[1, \\ \frac{\pi(\mathbf{r}')}{\pi(\mathbf{r}\_i)}
\\ z^{N-1} \right\]. \\

The \\z^{N-1}\\ factor (rather than \\z^N\\) arises because the stretch
is one-dimensional even though the parameter space is \\N\\-dimensional;
the Jacobian of the map \\\mathbf{r}\_i \mapsto
\mathbf{r}\_{\text{ref}} + z(\mathbf{r}\_i - \mathbf{r}\_{\text{ref}})\\
is \\z^N\\, but the proposal density contributes a factor of \\1/z\\
that partially cancels [\[1\]](#ref1).

### Directed proposals: the parabolic energy fit

The DQMC method augments QMC with energy guidance. After computing the
three QMC knots \\(t_i, t_j, t_k)\\ and evaluating the negative
log-posterior \\E = -\log\pi\\ at the corresponding positions, DQMC fits
a parabola \\E(t) = At^2 + Bt + C\\ through the three energy values by
solving the \\3 \times 3\\ interpolation system. When \\A \> 0\\ (the
parabola opens upward), the minimum is at \\t\_{\min} = -B/(2A)\\, and
the proposal knot is drawn from \\t' \sim \mathcal{N}(t\_{\min}, b^2)\\
rather than from the undirected distribution \\P_S\\. The parameter
\\b\\ controls how aggressively the proposal targets the energy minimum:
small \\b\\ concentrates proposals near \\t\_{\min}\\, while large \\b\\
diffuses them and recovers undirected-like behavior.

Because the forward and reverse proposal densities now differ, DQMC
requires an asymmetric Metropolis-Hastings correction. The acceptance
probability is

\\ \alpha = \min\left\[1, \\
\frac{\pi(\mathbf{R}(t'))}{\pi(\mathbf{r}\_i)} \\ \|w_i\|^N \\
\frac{P_S(t') \\ P_G(t_i \mid \text{new})}{P_S(t_i) \\ P_G(t' \mid
\text{old})} \right\], \\

where \\P_G(t \| \cdot)\\ is the Gaussian directed proposal density
centered on the parabolic minimum computed from the relevant set of
positions. The ratio \\P_S(t')/P_S(t_i)\\ corrects for the baseline knot
sampling, while \\P_G(t_i \| \text{new})/P_G(t' \| \text{old})\\
accounts for the asymmetry in the energy-directed step.

  

## The eight methods

### QMC (Quadratic Monte Carlo)

The basic QMC move selects two guide walkers \\j\\ and \\k\\ uniformly
at random from the complementary ensemble and constructs a second-order
Lagrange interpolant through the target walker \\i\\ and the two guides.
The target receives a random knot \\t_i \sim P_S(t)\\, the guides
receive fixed knots \\t_j = -1\\ and \\t_k = +1\\, and a proposal knot
\\t' \sim P_S(t)\\ generates weights \\(w_i, w_j, w_k)\\ that define the
proposal \\\mathbf{r}' = w_i \mathbf{r}\_i + w_j \mathbf{r}\_j + w_k
\mathbf{r}\_k\\. Acceptance uses the Jacobian \\\|w_i\|^N\\. This is the
foundation of the QMC family and serves as a general-purpose replacement
for AIES on curved posteriors. It requires a minimum ensemble of \\N_c
\geq 3\\ walkers and works well with `a` in the range 0.3–1.0 depending
on dimension.

#### Algorithm

1.  **Initialize** \\N_c \geq 3\\ walkers with positions \\\\z_1,
    \ldots, z\_{N_c}\\\\, each of dimension \\N\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  **For** each walker \\i = 1, \ldots, N_c\\:
        - Select two guide walkers \\j, k \neq i\\ uniformly at random
        - Draw target knot \\t_i \sim P_S(t)\\ and proposal knot \\t'
          \sim P_S(t)\\; fix guide knots \\t_j = -1\\, \\t_k = +1\\
        - Compute Lagrange weights \\w_i, w_j, w_k\\ from the three
          knots
        - Propose \\z^\* = w_i z_i + w_j z_j + w_k z_k\\
        - Compute \\\alpha = \min\\\bigl(1,\\ \|w_i\|^N \cdot \pi(z^\*)
          / \pi(z_i)\bigr)\\
        - With probability \\\alpha\\, set \\z_i \leftarrow z^\*\\
3.  **Return** samples from the primary walker

### QMC3 (third-order Monte Carlo)

QMC3 extends the interpolation to third order by selecting three guide
walkers with fixed knots at \\-1\\, \\0\\, and \\+1\\. The four-point
cubic interpolant can follow S-shaped or asymmetrically curved manifolds
that a quadratic cannot match. The cost is an additional guide walker
per proposal and slightly more extreme Lagrange weights for the same
knot range, so somewhat smaller `a` may be needed. Minimum ensemble size
is \\N_c \geq 4\\.

#### Algorithm

1.  **Initialize** \\N_c \geq 4\\ walkers with positions \\\\z_1,
    \ldots, z\_{N_c}\\\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  **For** each walker \\i\\:
        - Select three guide walkers \\j, k, l \neq i\\ at random; fix
          knots \\t_j = -1\\, \\t_k = 0\\, \\t_l = +1\\
        - Draw \\t_i \sim P_S(t)\\ and proposal knot \\t' \sim P_S(t)\\
        - Compute degree-3 Lagrange weights \\(w_i, w_j, w_k, w_l)\\
        - Propose \\z^\* = w_i z_i + w_j z_j + w_k z_k + w_l z_l\\
        - Accept with \\\alpha = \min\\\bigl(1,\\ \|w_i\|^N \cdot
          \pi(z^\*) / \pi(z_i)\bigr)\\
3.  **Return** samples from the primary walker

### QMCN (\\N\\th-order Monte Carlo)

QMCN generalizes to arbitrary polynomial order controlled by the
`nOrder` parameter. The method selects `nOrder` guide walkers placed at
equispaced knots on \\\[-1, +1\]\\ and constructs the degree-`nOrder`
Lagrange interpolant. This flexibility allows matching complex posterior
geometry, but diminishing returns set in beyond order 4 or 5 because the
Lagrange weights oscillate increasingly at high orders (Runge-like
behavior in the weight distribution). The minimum ensemble size is \\N_c
\geq n\_{\text{order}} + 1\\.

#### Algorithm

1.  **Initialize** \\N_c \geq n\_{\text{order}} + 1\\ walkers
2.  **For** \\t = 1, \ldots, T\\:
    1.  **For** each walker \\i\\:
        - Select \\n\_{\text{order}}\\ guide walkers at random; place at
          equispaced knots on \\\[-1, +1\]\\
        - Draw \\t_i \sim P_S(t)\\ and proposal knot \\t' \sim P_S(t)\\
        - Compute degree-\\n\_{\text{order}}\\ Lagrange weights
        - Propose \\z^\* = \sum_m w_m z_m\\
        - Accept with \\\alpha = \min\\\bigl(1,\\ \|w_i\|^N \cdot
          \pi(z^\*) / \pi(z_i)\bigr)\\
3.  **Return** samples from the primary walker

### DQMC (directed Quadratic Monte Carlo)

DQMC adds an energy-guided component to the standard QMC proposal. After
evaluating the target density at three knot positions, it fits a
parabola through the energy values and centers the proposal distribution
on the parabola’s minimum. The parameter `b` controls the width of this
guided Gaussian: small `b` produces aggressive moves toward the energy
minimum and works best on unimodal posteriors with strong curvature,
while large `b` relaxes toward undirected QMC. The asymmetric
Metropolis-Hastings correction accounts for the difference between
forward and reverse proposal densities. Minimum ensemble size is \\N_c
\geq 3\\.

#### Algorithm

1.  **Initialize** \\N_c \geq 3\\ walkers
2.  **For** \\t = 1, \ldots, T\\:
    1.  **For** each walker \\i\\:
        - Select two guide walkers \\j, k\\; draw \\t_i \sim P_S(t)\\;
          fix \\t_j = -1\\, \\t_k = +1\\
        - Evaluate energies \\E_i = -\log\pi(z_i)\\, \\E_j =
          -\log\pi(z_j)\\, \\E_k = -\log\pi(z_k)\\
        - Fit parabola \\E(t) = At^2 + Bt + C\\ through the three energy
          values
        - **If** \\A \> 0\\: draw \\t' \sim \mathcal{N}(-B/(2A),\\
          b^2)\\ (directed)
        - **Else**: draw \\t' \sim P_S(t)\\ (undirected fallback)
        - Compute Lagrange weights and propose \\z^\* = w_i z_i + w_j
          z_j + w_k z_k\\
        - Accept with asymmetric MH ratio: \\\alpha = \min\\\bigl(1,\\
          \|w_i\|^N \cdot \frac{\pi(z^\*)}{\pi(z_i)} \cdot \frac{P_S(t')
          \\ P_G(t_i \mid \text{new})}{P_S(t_i) \\ P_G(t' \mid
          \text{old})}\bigr)\\
3.  **Return** samples from the primary walker

### SQMC (Quadratic Simplex Monte Carlo)

SQMC replaces the individual guide walkers in QMC with weighted averages
of groups of walkers. Two groups of `nGuidePoints` walkers are selected,
and each group is collapsed to a randomly weighted centroid using
Dirichlet-like weights controlled by `weightRange`. The standard QMC
interpolation then uses these averaged positions instead of individual
walkers. This averaging suppresses individual walker noise and can
improve proposal quality when the ensemble is large, at the cost of
requiring more walkers: \\N_c \geq 2 \cdot n_G + 1\\ where \\n_G\\ is
`nGuidePoints`.

#### Algorithm

1.  **Initialize** \\N_c \geq 2 n_G + 1\\ walkers
2.  **For** \\t = 1, \ldots, T\\:
    1.  **For** each walker \\i\\:
        - Select two groups of \\n_G\\ walkers from the complementary
          ensemble
        - For each group, draw Dirichlet-like weights (controlled by
          `weightRange`) and compute weighted centroid \\\bar{z}\_j\\,
          \\\bar{z}\_k\\
        - Perform standard QMC interpolation using \\z_i\\,
          \\\bar{z}\_j\\, \\\bar{z}\_k\\ as the three interpolation
          points
        - Accept with \\\alpha = \min\\\bigl(1,\\ \|w_i\|^N \cdot
          \pi(z^\*) / \pi(z_i)\bigr)\\
3.  **Return** samples from the primary walker

### MAMC (modified affine Monte Carlo)

MAMC uses two guide walkers to define both the stretch direction and the
reference point, unlike standard AIES which stretches from a single
companion. The proposal is \\\mathbf{r}' = \mathbf{r}\_j +
z(\mathbf{r}\_i - \mathbf{r}\_k)\\ where \\z\\ is drawn from the stretch
distribution \\g(z) \propto 1/\sqrt{z}\\ on \\\[1/a, a\]\\. By using two
guides to define the direction, MAMC can capture correlations between
guide positions that AIES misses. The acceptance correction is
\\z^{N-1}\\ as in standard AIES. The parameter `a` must be greater than
1 (default 2.0); minimum ensemble size is \\N_c \geq 3\\.

#### Algorithm

1.  **Initialize** \\N_c \geq 3\\ walkers
2.  **For** \\t = 1, \ldots, T\\:
    1.  **For** each walker \\i\\:
        - Select two guide walkers \\j, k \neq i\\ from the
          complementary ensemble
        - Draw stretch factor \\z \sim g(z) \propto 1/\sqrt{z}\\ on
          \\\[1/a,\\ a\]\\
        - Propose \\z^\* = z_j + z \cdot (z_i - z_k)\\
        - Accept with \\\alpha = \min\\\bigl(1,\\ z^{N-1} \cdot
          \pi(z^\*) / \pi(z_i)\bigr)\\
3.  **Return** samples from the primary walker

### SAMC (affine Simplex Monte Carlo)

SAMC combines the AIES stretch move with simplex averaging: instead of
stretching from a single companion, it stretches from the weighted
centroid of `nGuidePoints` walkers. The proposal is \\\mathbf{r}' =
\bar{\mathbf{r}} + z(\mathbf{r}\_i - \bar{\mathbf{r}})\\ where
\\\bar{\mathbf{r}}\\ is the guide centroid and \\z\\ follows the stretch
distribution. The centroid provides a more stable reference point than a
single walker, which can reduce proposal variance. The acceptance
correction is \\z^{N-1}\\; `a` must be greater than 1 (default 2.0);
minimum ensemble size is \\N_c \geq n_G + 1\\.

#### Algorithm

1.  **Initialize** \\N_c \geq n_G + 1\\ walkers
2.  **For** \\t = 1, \ldots, T\\:
    1.  **For** each walker \\i\\:
        - Select \\n_G\\ guide walkers from the complementary ensemble
        - Draw Dirichlet-like weights (controlled by `weightRange`) and
          compute centroid \\\bar{z}\\
        - Draw stretch factor \\z \sim g(z) \propto 1/\sqrt{z}\\ on
          \\\[1/a,\\ a\]\\
        - Propose \\z^\* = \bar{z} + z \cdot (z_i - \bar{z})\\
        - Accept with \\\alpha = \min\\\bigl(1,\\ z^{N-1} \cdot
          \pi(z^\*) / \pi(z_i)\bigr)\\
3.  **Return** samples from the primary walker

### WMC (Walk Monte Carlo)

WMC constructs a symmetric random walk perturbation from the ensemble
covariance structure. The proposal is \\\mathbf{r}' = \mathbf{r}\_i +
\sum\_{g=1}^{n_G} a \\ Z_g \\ (\mathbf{r}\_g - \bar{\mathbf{r}})\\ where
\\Z_g \sim \mathcal{N}(0,1)\\, \\\mathbf{r}\_g\\ are guide walker
positions, and \\\bar{\mathbf{r}}\\ is their mean. Because the Gaussian
perturbation is symmetric around \\\mathbf{r}\_i\\, no Jacobian
correction is needed; the acceptance reduces to the standard Metropolis
ratio \\\min\[1, \pi(\mathbf{r}')/\pi(\mathbf{r}\_i)\]\\. WMC is the
simplest method in the family and serves as a robust baseline when the
other methods’ Jacobian penalties become problematic. Minimum ensemble
size is \\N_c \geq n_G + 1\\.

#### Algorithm

1.  **Initialize** \\N_c \geq n_G + 1\\ walkers
2.  **For** \\t = 1, \ldots, T\\:
    1.  **For** each walker \\i\\:
        - Select \\n_G\\ guide walkers \\\\z\_{g_1}, \ldots,
          z\_{g\_{n_G}}\\\\; compute their mean \\\bar{z}\\
        - Draw \\Z_g \sim \mathcal{N}(0, 1)\\ for \\g = 1, \ldots, n_G\\
        - Propose \\z^\* = z_i + \sum\_{g=1}^{n_G} a \\ Z_g \\
          (z\_{g_g} - \bar{z})\\
        - Accept with \\\alpha = \min\\\bigl(1,\\ \pi(z^\*) /
          \pi(z_i)\bigr)\\ (no Jacobian correction)
3.  **Return** samples from the primary walker

  

## Benchmark problems

The following four distributions serve as test problems for comparing
the QMC methods. Each has known analytical properties that allow
quantitative assessment of posterior recovery.

### Multivariate normal (5D)

A 5-dimensional standard normal distribution provides a basic sanity
check. The target is \\\pi(\mathbf{x}) = \mathcal{N}(\mathbf{0},
\mathbf{I}\_5)\\, so \\\log \pi(\mathbf{x}) =
-\frac{1}{2}\\\mathbf{x}\\^2\\ up to a normalizing constant. The known
posterior mean is \\\mathbf{0}\\ and the marginal variance is 1 for
every coordinate.

``` r

NormalModel <- function(parm, Data) {
  LP <- -0.5 * sum(parm^2)
  list(LP = LP, Dev = -2 * LP, Monitor = LP,
       yhat = rep(0, Data$N), parm = parm)
}

NormalData <- list(
  mon.names = "LP",
  parm.names = paste0("x", 1:5),
  N = 1
)
```

### Rosenbrock banana (2D)

The Rosenbrock banana is a bivariate distribution with strong nonlinear
correlation:

\\ \log \pi(x_1, x_2) = -\frac{1}{2}\left\[\frac{x_1^2}{\sigma_1^2} +
\frac{(x_2 - b x_1^2 + c)^2}{\sigma_2^2}\right\], \\

with \\b = 0.1\\, \\c = 10\\, \\\sigma_1 = 10\\, \\\sigma_2 = 1\\. The
marginal \\x_1 \sim \mathcal{N}(0, 100)\\ is straightforward, but the
conditional \\x_2 \mid x_1 \sim \mathcal{N}(b x_1^2 - c, 1)\\ creates
the characteristic curved ridge. The analytical moments are
\\\mathbb{E}\[x_1\] = 0\\, \\\mathbb{E}\[x_2\] = b\sigma_1^2 - c = 0\\,
\\\text{Var}(x_1) = 100\\, and \\\text{Var}(x_2) = \sigma_2^2 +
2b^2\sigma_1^4 = 201\\. This distribution is the natural testbed for QMC
methods because their quadratic proposals can track the parabolic ridge.

``` r

RosenbrockModel <- function(parm, Data) {
  x1 <- parm[1]; x2 <- parm[2]
  LP <- -0.5 * (x1^2 / 100 + (x2 - 0.1 * x1^2 + 10)^2)
  list(LP = LP, Dev = -2 * LP, Monitor = LP,
       yhat = rep(0, Data$N), parm = parm)
}

RosenbrockData <- list(
  mon.names = "LP",
  parm.names = c("x1", "x2"),
  N = 1
)
```

### Ring distribution (2D)

The ring distribution concentrates mass on an annulus of radius \\R\\:

\\ \log \pi(x_1, x_2) = -\frac{(\rho - R)^{2m}}{2T}, \qquad \rho =
\sqrt{x_1^2 + x_2^2}, \\

with \\R = 5\\, \\m = 2\\, \\T = 0.5\\. The quartic potential \\(\rho -
5)^4\\ creates a narrow ring (density drops to \\e^{-1}\\ of maximum at
\\\rho \approx 5 \pm 1\\). The distribution is rotationally symmetric,
so \\\mathbb{E}\[x_1\] = \mathbb{E}\[x_2\] = 0\\ and
\\\mathbb{E}\[\rho\] \approx R = 5\\. The ring geometry is challenging
because the walker population must spread around the full \\2\pi\\
angular range; methods that propose moves along the ring (rather than
radially) will mix faster.

``` r

RingModel <- function(parm, Data) {
  rho <- sqrt(parm[1]^2 + parm[2]^2)
  LP <- -((rho - 5)^4) / (2 * 0.5)
  list(LP = LP, Dev = -2 * LP, Monitor = LP,
       yhat = rep(0, Data$N), parm = parm)
}

RingData <- list(
  mon.names = "LP",
  parm.names = c("x1", "x2"),
  N = 1
)
```

### High-dimensional banana (8D)

This distribution couples four pairs of coordinates through the same
nonlinear banana transformation, extending the Rosenbrock geometry to 8
dimensions:

\\ \log \pi(\mathbf{x}) = -\frac{1}{2} \sum\_{k=1}^{4}
\left\[\frac{x\_{2k-1}^2}{100} + (x\_{2k} + 0.1 \\ x\_{2k-1}^2 -
10)^2\right\]. \\

The odd-indexed coordinates have marginal distributions \\x\_{2k-1} \sim
\mathcal{N}(0, 100)\\, and the even-indexed coordinates are
conditionally normal given their paired odd coordinate. The
8-dimensional setting tests how well the Jacobian penalty \\\|w_i\|^N\\
or \\z^{N-1}\\ scales with dimension.

``` r

BananaModel <- function(parm, Data) {
  d <- length(parm)
  LP <- 0
  for (i in seq(1, d - 1, by = 2)) {
    x1 <- parm[i]; x2 <- parm[i + 1]
    y2 <- x2 + 0.1 * x1^2 - 10
    LP <- LP - 0.5 * (x1^2 / 100 + y2^2)
  }
  list(LP = LP, Dev = -2 * LP, Monitor = LP,
       yhat = rep(0, Data$N), parm = parm)
}

BananaData <- list(
  mon.names = "LP",
  parm.names = paste0("x", 1:8),
  N = 1
)
```

  

## Comparative benchmarks

### Helper functions

The `make_specs` function constructs the method-appropriate `Specs`
list, and `run_one` runs a single benchmark and collects the fitted
`demonoid` object alongside the post-burn-in chain and diagnostic
summaries computed with lucifer’s modern diagnostic functions (bulk-ESS,
tail-ESS, split \\\hat{R}\\).

``` r

make_specs <- function(method, Nc = 20, a_qmc = 0.6, a_affine = 2.0) {
  switch(method,
    QMC  = list(Nc = Nc, Z = NULL, a = a_qmc, Gaussian = FALSE),
    QMC3 = list(Nc = Nc, Z = NULL, a = a_qmc, Gaussian = FALSE),
    QMCN = list(Nc = Nc, Z = NULL, a = a_qmc, Gaussian = FALSE, nOrder = 4L),
    DQMC = list(Nc = Nc, Z = NULL, a = a_qmc, Gaussian = FALSE, b = 1.0),
    SQMC = list(Nc = Nc, Z = NULL, a = a_qmc, Gaussian = FALSE,
                nGuidePoints = 3L, weightRange = 1.0),
    MAMC = list(Nc = Nc, Z = NULL, a = a_affine),
    SAMC = list(Nc = Nc, Z = NULL, a = a_affine,
                nGuidePoints = 3L, weightRange = 1.0),
    WMC  = list(Nc = Nc, Z = NULL, a = a_qmc, nGuidePoints = 3L),
    AIES = list(Nc = Nc, Z = NULL, beta = a_affine,
                CPUs = 1, Packages = NULL, Dyn.libs = NULL)
  )
}

run_one <- function(method, Model, Data, iv, n_iter = 20000,
                    Nc = 20, a_qmc = 0.6, a_affine = 2.0, seed = 42) {
  burn <- floor(n_iter * 0.25)
  specs <- make_specs(method, Nc, a_qmc, a_affine)
  set.seed(seed)
  invisible(capture.output({
    fit <- lucifer(Model, Data, Initial.Values = iv,
      Iterations = n_iter, Status = n_iter + 1, Thinning = 1,
      Algorithm = method, Chains = 1, Specs = specs)
  }))
  chain <- fit$Posterior1[(burn + 1):n_iter, , drop = FALSE]
  list(
    fit = fit,
    acceptance = fit$Acceptance.Rate,
    chain = chain,
    post_mean = colMeans(chain),
    post_sd = apply(chain, 2, sd),
    ess_bulk = ESS.bulk(chain),
    ess_tail = ESS.tail(chain),
    rhat = Rhat(chain)
  )
}

all_methods <- c("QMC", "QMC3", "QMCN", "DQMC", "SQMC",
                 "MAMC", "SAMC", "WMC", "AIES")
```

### Multivariate normal (5D)

``` r

results_normal <- lapply(all_methods, function(m) {
  run_one(m, NormalModel, NormalData, iv = rep(0, 5),
          n_iter = 30000, a_qmc = 0.6, a_affine = 2.0)
})
names(results_normal) <- all_methods
```

The true posterior is \\\mathcal{N}(\mathbf{0}, \mathbf{I}\_5)\\, so all
marginal means are zero and all marginal standard deviations are 1. The
table reports acceptance rate, posterior mean and SD for the first
coordinate (representative, since the target is isotropic), mean
bulk-ESS across coordinates, and the corresponding ESS per post-burn-in
iteration.

``` r

tab_normal <- data.frame(
  Method = all_methods,
  `Accept.` = sapply(results_normal, function(r) round(r$acceptance, 3)),
  `E[x1]` = sapply(results_normal, function(r) round(r$post_mean[1], 3)),
  `SD(x1)` = sapply(results_normal, function(r) round(r$post_sd[1], 3)),
  `Bulk ESS` = sapply(results_normal, function(r) round(mean(r$ess_bulk), 0)),
  `ESS/iter` = sapply(results_normal, function(r)
    round(mean(r$ess_bulk) / nrow(r$chain), 3)),
  check.names = FALSE, row.names = NULL
)
knitr::kable(tab_normal,
  caption = "5D standard normal (true: E[x]=0, SD=1)",
  align = "lrrrrr")
```

On this isotropic target all methods achieve reasonable acceptance rates
and recover the true moments accurately. The ESS differences reflect how
efficiently each proposal explores the 5-dimensional space. WMC, which
carries no Jacobian correction, tends to produce the highest ESS per
iteration on symmetric targets. The affine methods (MAMC, SAMC, AIES)
benefit from the mild \\z^{N-1}\\ penalty in low dimensions, while
QMC-type methods pay the \\\|w_i\|^N\\ cost without gaining the
curvature advantage that makes them shine on non-Gaussian targets.

### Rosenbrock banana (2D)

``` r

results_rosen <- lapply(all_methods, function(m) {
  run_one(m, RosenbrockModel, RosenbrockData, iv = c(0, 0),
          n_iter = 50000, a_qmc = 0.8, a_affine = 2.0)
})
names(results_rosen) <- all_methods
```

The analytical moments are \\\mathbb{E}\[x_1\] = 0\\,
\\\mathbb{E}\[x_2\] = 0\\, \\\text{SD}(x_1) = 10\\, \\\text{SD}(x_2) =
\sqrt{201} \approx 14.18\\. The table includes the true values as a
reference row.

``` r

true_row_rosen <- data.frame(
  Method = "True", `Accept.` = NA,
  `E[x1]` = 0, `E[x2]` = 0,
  `SD(x1)` = 10.00, `SD(x2)` = round(sqrt(201), 2),
  `Bulk ESS` = NA, check.names = FALSE)
tab_rosen <- data.frame(
  Method = all_methods,
  `Accept.` = sapply(results_rosen, function(r) round(r$acceptance, 3)),
  `E[x1]` = sapply(results_rosen, function(r) round(r$post_mean[1], 2)),
  `E[x2]` = sapply(results_rosen, function(r) round(r$post_mean[2], 2)),
  `SD(x1)` = sapply(results_rosen, function(r) round(r$post_sd[1], 2)),
  `SD(x2)` = sapply(results_rosen, function(r) round(r$post_sd[2], 2)),
  `Bulk ESS` = sapply(results_rosen, function(r) round(mean(r$ess_bulk), 0)),
  check.names = FALSE
)
tab_rosen <- rbind(true_row_rosen, tab_rosen)
knitr::kable(tab_rosen, row.names = FALSE,
  caption = "2D Rosenbrock banana: posterior recovery and efficiency",
  align = "lrrrrrr")
```

The Rosenbrock banana is where the QMC family is designed to shine.
Quadratic interpolation naturally tracks the parabolic ridge \\x_2
\approx 0.1 x_1^2 - 10\\, producing proposals that stay in the
high-density region rather than cutting across it. Methods that follow
the curvature (QMC, QMC3, DQMC) should show higher ESS than AIES, whose
linear stretch moves overshoot the banana ridge and waste proposals.
AIES is expected to underperform here because its straight-line
interpolation between two walkers crosses perpendicular to the narrow
curved manifold, leading to high rejection rates and persistent
autocorrelation. This is precisely the geometric limitation that
motivated the QMC family.

### Ring distribution (2D)

``` r

results_ring <- lapply(all_methods, function(m) {
  run_one(m, RingModel, RingData, iv = c(5, 0),
          n_iter = 50000, a_qmc = 0.8, a_affine = 2.0)
})
names(results_ring) <- all_methods
```

The ring is centered at the origin with radius \\R = 5\\. By rotational
symmetry, \\\mathbb{E}\[x_1\] = \mathbb{E}\[x_2\] = 0\\ and
\\\mathbb{E}\[\rho\] \approx 5\\.

``` r

true_row_ring <- data.frame(
  Method = "True", `Accept.` = NA,
  `E[x1]` = 0, `E[x2]` = 0,
  `E[rho]` = 5.00,
  `Bulk ESS` = NA, `ESS/iter` = NA,
  check.names = FALSE)
tab_ring <- data.frame(
  Method = all_methods,
  `Accept.` = sapply(results_ring, function(r) round(r$acceptance, 3)),
  `E[x1]` = sapply(results_ring, function(r) round(r$post_mean[1], 2)),
  `E[x2]` = sapply(results_ring, function(r) round(r$post_mean[2], 2)),
  `E[rho]` = sapply(results_ring, function(r) {
    chain <- r$chain
    round(mean(sqrt(chain[, 1]^2 + chain[, 2]^2)), 2)
  }),
  `Bulk ESS` = sapply(results_ring, function(r) round(mean(r$ess_bulk), 0)),
  `ESS/iter` = sapply(results_ring, function(r)
    round(mean(r$ess_bulk) / nrow(r$chain), 3)),
  check.names = FALSE
)
tab_ring <- rbind(true_row_ring, tab_ring)
knitr::kable(tab_ring, row.names = FALSE,
  caption = "2D ring (R=5): posterior recovery and efficiency",
  align = "lrrrrrr")
```

The ring distribution tests angular exploration rather than curvature
tracking. All walkers start near \\(5, 0)\\ and must spread around the
full circle during burn-in. None of the QMC-type methods has a geometric
advantage here because the ring’s curvature is circular rather than
parabolic, meaning that polynomial proposals along the chord are no
better than linear ones for traversing the angular direction. Methods
that produce diverse proposal directions (WMC, SAMC with multiple
guides) may explore the ring more efficiently than single-direction
methods, and AIES stretch moves can perform well on the ring because the
radial concentration of mass is simple even though the angular extent is
challenging.

### High-dimensional banana (8D)

``` r

results_banana <- lapply(all_methods, function(m) {
  run_one(m, BananaModel, BananaData, iv = rep(0, 8),
          n_iter = 30000, a_qmc = 0.3, a_affine = 2.0)
})
names(results_banana) <- all_methods
```

The odd-indexed coordinates have true marginal mean 0 and SD 10; the
even-indexed coordinates have conditional mean \\10 - 0.1 x\_{2k-1}^2\\
given the odd partner, with marginal mean 0 and SD \\\sqrt{201} \approx
14.18\\.

``` r

true_row_banana <- data.frame(
  Method = "True", `Accept.` = NA,
  `E[x]` = 0,
  `SD(odd)` = 10.00, `SD(even)` = round(sqrt(201), 2),
  `Bulk ESS` = NA, `ESS/iter` = NA,
  check.names = FALSE)
tab_banana <- data.frame(
  Method = all_methods,
  `Accept.` = sapply(results_banana, function(r) round(r$acceptance, 3)),
  `E[x]` = sapply(results_banana, function(r)
    round(mean(r$post_mean), 2)),
  `SD(odd)` = sapply(results_banana, function(r)
    round(mean(r$post_sd[seq(1, 7, by = 2)]), 2)),
  `SD(even)` = sapply(results_banana, function(r)
    round(mean(r$post_sd[seq(2, 8, by = 2)]), 2)),
  `Bulk ESS` = sapply(results_banana, function(r) round(mean(r$ess_bulk), 0)),
  `ESS/iter` = sapply(results_banana, function(r)
    round(mean(r$ess_bulk) / nrow(r$chain), 3)),
  check.names = FALSE, row.names = NULL
)
tab_banana <- rbind(true_row_banana, tab_banana)
knitr::kable(tab_banana, row.names = FALSE,
  caption = "8D banana: posterior recovery and efficiency (a_qmc=0.3)",
  align = "lrrrrrr")
```

In 8 dimensions, the Jacobian penalty \\\|w_i\|^8\\ is substantially
harsher than \\\|w_i\|^2\\, requiring smaller `a` for QMC-type methods
(`a = 0.3` here vs `a = 0.8` in 2D). The `SD(odd)` and `SD(even)`
columns report the mean estimated standard deviation across the four
odd-indexed (\\x_1, x_3, x_5, x_7\\, true SD = 10) and four even-indexed
(\\x_2, x_4, x_6, x_8\\, true SD \\\approx\\ 14.18) coordinates
respectively. Because the even-indexed coordinates involve the nonlinear
banana transformation \\x\_{2k} \mid x\_{2k-1} \sim \mathcal{N}(10 - 0.1
x\_{2k-1}^2, 1)\\, the SD(even) column is particularly sensitive to
underexploration of the banana tails. The smaller QMC step sizes reduce
the curvature advantage, and the relative performance of different
methods becomes more tightly clustered, with WMC often emerging as the
most robust choice in higher dimensions because it avoids the
exponentially growing Jacobian penalty entirely.

### Summary comparison

The following table condenses the benchmark results into acceptance rate
and ESS rank (1 = best) for each problem, along with the average rank
across all four benchmarks.

``` r

get_rank <- function(results) {
  ess_vals <- sapply(results, function(r) mean(r$ess_bulk))
  rank(-ess_vals, ties.method = "min")
}
ranks <- data.frame(
  Method = all_methods,
  `Normal AR` = sapply(results_normal, function(r) round(r$acceptance, 2)),
  `Normal rank` = get_rank(results_normal),
  `Rosen. AR` = sapply(results_rosen, function(r) round(r$acceptance, 2)),
  `Rosen. rank` = get_rank(results_rosen),
  `Ring AR` = sapply(results_ring, function(r) round(r$acceptance, 2)),
  `Ring rank` = get_rank(results_ring),
  `Banana AR` = sapply(results_banana, function(r) round(r$acceptance, 2)),
  `Banana rank` = get_rank(results_banana),
  check.names = FALSE, row.names = NULL
)
ranks$`Avg rank` <- round(rowMeans(ranks[, c(3, 5, 7, 9)]), 1)
knitr::kable(ranks, row.names = FALSE,
  caption = "Summary: acceptance rates and ESS ranks across benchmarks",
  align = "lrrrrrrrr")
```

No single method dominates across all problems. QMC-type methods tend to
rank well on curved posteriors (Rosenbrock, banana) where their
polynomial proposals match the geometry, while the simpler methods (WMC,
AIES) are competitive on isotropic targets (normal) and angular
exploration (ring). Several patterns are worth noting. First, acceptance
rate does not predict ESS rank: a method can accept many proposals that
are highly correlated (high acceptance, low ESS), or fewer proposals
that are nearly independent (moderate acceptance, high ESS). Second, the
QMC family’s advantage is geometry-dependent; on the isotropic normal or
the ring, the curvature-tracking capability is irrelevant, and the
Jacobian cost becomes a net liability. Third, WMC is the most robust
all-around performer because its lack of Jacobian correction means it
never pays an exponential penalty for dimension, making it a safe
default when the posterior geometry is unknown. The overall pattern
supports Militzer’s recommendation [\[4\]](#ref4) to try QMC as the
first alternative when AIES underperforms on curved posteriors, and to
fall back to WMC in high dimensions.

  

## Diagnostic analysis

The QMC ensemble methods produce standard `demonoid` objects, so all of
lucifer’s diagnostic infrastructure applies directly. This section
demonstrates the available tools and provides a comprehensive comparison
of convergence quality across methods.

### Built-in diagnostics

Every `demonoid` object returned by
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
supports [`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html), and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods. The
[`summary()`](https://rdrr.io/r/base/summary.html) method computes
acceptance rate, stationarity assessment (via Hellinger-distance batch
diagnostics), DIC, split \\\hat{R}\\, bulk-ESS, tail-ESS, and MCSE for
all monitored parameters. For a representative example, the following
shows the summary output for QMC on the Rosenbrock banana.

``` r

# Available for any fitted demonoid object:
summary(results_rosen[["QMC"]]$fit)

# Built-in trace/density/ACF plots (3-column grid per parameter):
plot(results_rosen[["QMC"]]$fit, Data = RosenbrockData)
```

The
[`plot.demonoid()`](https://robustecologies.github.io/lucifer/reference/plot.demonoid.md)
method produces a three-column panel for each parameter: trace plot,
kernel density estimate, and autocorrelation function, identical to the
custom plots shown below but applied automatically to all parameters
including monitored quantities. The
[`Consort()`](https://robustecologies.github.io/lucifer/reference/Consort.md)
function provides a text-based diagnostic advisor that evaluates five
convergence criteria (non-adaptive algorithm, acceptance rate, MCSE,
ESS, and stationarity) and recommends next steps.

### Convergence: split \\\hat{R}\\ and Heidelberger diagnostic

Split \\\hat{R}\\ [\[8\]](#ref8) measures within-chain convergence by
splitting the chain in half and comparing the two halves; values below
1.01 indicate excellent convergence, while values above 1.05 warrant
concern. The Heidelberger diagnostic applies a Cramer-von Mises
stationarity test and a halfwidth criterion to each marginal posterior.

``` r

# Wrapper to run Heidelberger on post-burn-in chain (Posterior2 may be too
# short because the automatic burn-in recommendation can consume nearly all
# samples on ensemble methods with short initial transients)
run_heidelberger <- function(chain) {
  obj <- list(Posterior1 = chain, Posterior2 = NA)
  class(obj) <- "demonoid"
  Heidelberger.Diagnostic(obj)
}

heid_pass <- function(results, method) {
  hd <- tryCatch(run_heidelberger(results[[method]]$chain),
                 error = function(e) NULL)
  if (is.null(hd)) return("n/a")
  paste0(sum(hd[, "stest"] == 1), "/", nrow(hd))
}

conv_all <- data.frame(
  Method = all_methods,
  `Rosen. Rhat` = sapply(results_rosen, function(r)
    round(max(r$rhat, na.rm = TRUE), 4)),
  `Rosen. Heid.` = sapply(all_methods, function(m) heid_pass(results_rosen, m)),
  `Ring Rhat` = sapply(results_ring, function(r)
    round(max(r$rhat, na.rm = TRUE), 4)),
  `Ring Heid.` = sapply(all_methods, function(m) heid_pass(results_ring, m)),
  `Banana Rhat` = sapply(results_banana, function(r)
    round(max(r$rhat, na.rm = TRUE), 4)),
  `Banana Heid.` = sapply(all_methods, function(m) heid_pass(results_banana, m)),
  check.names = FALSE, row.names = NULL
)
knitr::kable(conv_all, row.names = FALSE,
  caption = "Convergence diagnostics across benchmarks (Heidelberger: parameters passing stationarity test)",
  align = "lcccccc")
```

Methods with \\\hat{R}\\ close to 1.00 and all parameters passing the
Heidelberger stationarity test have achieved adequate convergence within
the 20,000-iteration budget. Methods with elevated \\\hat{R}\\ or failed
Heidelberger tests need longer runs or better tuning. AIES is expected
to show weaker convergence on the Rosenbrock banana because its linear
proposals produce high autocorrelation on curved manifolds, requiring
substantially more iterations to achieve stationarity.

### Effective sample size

Bulk-ESS measures mixing quality for central posterior summaries (means,
medians), while tail-ESS measures mixing for extreme quantiles (credible
interval endpoints). Both should exceed 400 for reliable inference
[\[8\]](#ref8). The table below reports the mean across parameters for
each method-benchmark pair.

``` r

ess_tab <- data.frame(
  Method = all_methods,
  `Rosen. Bulk` = sapply(results_rosen, function(r) round(mean(r$ess_bulk), 0)),
  `Rosen. Tail` = sapply(results_rosen, function(r) round(mean(r$ess_tail), 0)),
  `Ring Bulk` = sapply(results_ring, function(r) round(mean(r$ess_bulk), 0)),
  `Ring Tail` = sapply(results_ring, function(r) round(mean(r$ess_tail), 0)),
  `Banana Bulk` = sapply(results_banana, function(r) round(mean(r$ess_bulk), 0)),
  `Banana Tail` = sapply(results_banana, function(r) round(mean(r$ess_tail), 0)),
  check.names = FALSE, row.names = NULL
)
knitr::kable(ess_tab, row.names = FALSE,
  caption = "Bulk-ESS and tail-ESS across benchmarks (mean over parameters)",
  align = "lrrrrrr")
```

Tail-ESS is typically lower than bulk-ESS because the tails of the
posterior are explored less frequently. Methods with tail-ESS well below
400 may produce unreliable credible intervals even when their bulk-ESS
is adequate. On the banana targets, the tail-ESS penalty is most severe
for AIES because the linear proposals fail to reach the extreme arms of
the banana distribution, while QMC methods with their curvature-tracking
proposals can follow the ridge into the tails.

### Monte Carlo standard error

The MCSE quantifies the precision of the posterior mean estimate:
smaller MCSE means the estimate is closer to the true posterior mean
given infinite samples. The MCSE should be small relative to the
posterior standard deviation; a common criterion is MCSE \< 5% of SD
[\[9\]](#ref9).

``` r

mcse_rosen <- data.frame(
  Method = all_methods,
  `MCSE(x1)` = sapply(all_methods, function(m) {
    round(MCSE(results_rosen[[m]]$chain[, 1]), 3)
  }),
  `SD(x1)` = sapply(results_rosen, function(r) round(r$post_sd[1], 2)),
  `MCSE/SD` = sapply(all_methods, function(m) {
    mc <- MCSE(results_rosen[[m]]$chain[, 1])
    sd_val <- sd(results_rosen[[m]]$chain[, 1])
    paste0(round(100 * mc / sd_val, 1), "%")
  }),
  `MCSE(x2)` = sapply(all_methods, function(m) {
    round(MCSE(results_rosen[[m]]$chain[, 2]), 3)
  }),
  `SD(x2)` = sapply(results_rosen, function(r) round(r$post_sd[2], 2)),
  check.names = FALSE, row.names = NULL
)
knitr::kable(mcse_rosen, row.names = FALSE,
  caption = "MCSE for Rosenbrock banana parameters (MCSE/SD < 5% is adequate)",
  align = "lrrcrr")
```

Methods with high autocorrelation produce larger MCSE for the same chain
length, meaning longer runs are needed to achieve the same precision.
This is the practical cost of poor mixing: not incorrect estimates per
se, but imprecise ones that require more computation to narrow.

### Trace plots: Rosenbrock \\x_1\\

``` r

diag_methods <- c("QMC", "DQMC", "WMC", "AIES")
method_colors <- c(QMC = "#440154", DQMC = "#31688e",
                   WMC = "#35b779", AIES = "#fde725")
```

``` r

traces_rosen <- do.call(rbind, lapply(diag_methods, function(m) {
  chain <- results_rosen[[m]]$chain
  data.frame(
    iteration = seq_len(nrow(chain)),
    x1 = chain[, 1],
    method = m
  )
}))
traces_rosen$method <- factor(traces_rosen$method, levels = diag_methods)

ggplot(traces_rosen, aes(x = iteration, y = x1, color = method)) +
  geom_line(alpha = 0.5, linewidth = 0.3) +
  facet_wrap(~method, ncol = 2) +
  scale_color_manual(values = method_colors) +
  labs(x = "Iteration (post burn-in)", y = expression(x[1]),
       title = "Rosenbrock banana: trace plots") +
  theme_bw() +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))
```

Good mixing appears as rapid, uncorrelated oscillations filling the
marginal range \\x_1 \in \[-30, 30\]\\ (approximately \\\pm 3\sigma\\
for \\\sigma = 10\\). Persistent trends or slow drifts indicate high
autocorrelation and poor exploration of the curved manifold. AIES
typically shows the slowest mixing on this target because its linear
stretch moves cut perpendicular to the banana ridge, causing the walker
to take small effective steps along the ridge and producing long
autocorrelation tails. QMC and DQMC, by contrast, propose along the
curved manifold and achieve faster traversal of the full \\x_1\\ range.

### Posterior density: Rosenbrock \\x_1\\

``` r

# Compute kernel density estimates for each method
dens_list <- lapply(diag_methods, function(m) {
  d <- density(results_rosen[[m]]$chain[, 1], n = 512)
  data.frame(x1 = d$x, density = d$y, method = m)
})
dens_rosen <- do.call(rbind, dens_list)
dens_rosen$method <- factor(dens_rosen$method, levels = diag_methods)

# Analytical reference on the same grid
x_grid <- seq(-40, 40, length.out = 512)
ref_df <- data.frame(x1 = x_grid, density = dnorm(x_grid, 0, 10))

ggplot(dens_rosen, aes(x = x1, y = density, color = method)) +
  geom_line(linewidth = 0.7) +
  geom_line(data = ref_df, aes(x = x1, y = density),
            inherit.aes = FALSE, color = "black",
            linewidth = 1, linetype = "dashed") +
  scale_y_log10() +
  scale_color_manual(values = method_colors) +
  labs(x = expression(x[1]), y = expression(paste("Density (", log[10], " scale)")),
       title = expression(paste("Rosenbrock: marginal density of ", x[1],
                                " vs analytical ", N(0, 100)))) +
  theme_bw() +
  theme(legend.position = "bottom")
```

The dashed black line shows the true marginal \\\mathcal{N}(0, 100)\\.
The log-scale y-axis reveals differences in the tails that are invisible
on a linear scale: methods that underexplore the banana arms show
densities that fall off too steeply relative to the Gaussian reference,
while methods that get stuck in one region show localized spikes orders
of magnitude above the true density. AIES typically shows a very narrow,
high spike because the walker remains confined to a small region of the
banana ridge, producing a degenerate density estimate far from the true
Gaussian shape.

### Two-dimensional posterior: Rosenbrock

``` r

scatter_rosen <- do.call(rbind, lapply(diag_methods, function(m) {
  chain <- results_rosen[[m]]$chain
  data.frame(x1 = chain[, 1], x2 = chain[, 2], method = m)
}))
scatter_rosen$method <- factor(scatter_rosen$method, levels = diag_methods)

ggplot(scatter_rosen, aes(x = x1, y = x2)) +
  geom_point(alpha = 0.1, size = 0.3, color = "steelblue") +
  stat_density_2d(color = "darkblue", linewidth = 0.4, bins = 8) +
  facet_wrap(~method, ncol = 2) +
  labs(x = expression(x[1]), y = expression(x[2]),
       title = "Rosenbrock banana: 2D posterior samples") +
  theme_bw() +
  theme(strip.text = element_text(face = "bold"))
```

The posterior concentrates along the parabolic ridge \\x_2 \approx 0.1
x_1^2 - 10\\. Methods that adequately cover this ridge show the
characteristic banana shape spanning the full range \\x_1 \in \[-30,
+30\]\\, with smooth contour lines following the parabola. AIES may show
truncated or asymmetric coverage because its linear proposals cross the
banana ridge rather than following it; the walker effectively performs a
constrained random walk along a narrow channel, and 15,000 post-burn-in
iterations may be insufficient to traverse the full extent. QMC and DQMC
panels should show broader, more symmetric coverage because their
quadratic proposals track the curvature naturally.

### Two-dimensional posterior: ring

``` r

scatter_ring <- do.call(rbind, lapply(diag_methods, function(m) {
  chain <- results_ring[[m]]$chain
  data.frame(x1 = chain[, 1], x2 = chain[, 2], method = m)
}))
scatter_ring$method <- factor(scatter_ring$method, levels = diag_methods)

ggplot(scatter_ring, aes(x = x1, y = x2)) +
  geom_point(alpha = 0.1, size = 0.3, color = "firebrick") +
  stat_density_2d(color = "darkred", linewidth = 0.4, bins = 8) +
  facet_wrap(~method, ncol = 2) +
  coord_equal(xlim = c(-8, 8), ylim = c(-8, 8)) +
  labs(x = expression(x[1]), y = expression(x[2]),
       title = "Ring distribution: 2D posterior samples") +
  theme_bw() +
  theme(strip.text = element_text(face = "bold"))
```

Good exploration of the ring appears as a uniform distribution of
samples around the circle at radius 5. Gaps or clustering at specific
angles indicate that the sampler failed to explore the full angular
range during the run. All walkers start near \\(5, 0)\\, so coverage of
the left semicircle (\\x_1 \< 0\\) is a particularly informative
diagnostic for mixing speed. The ring geometry does not favor polynomial
proposals over linear ones: the curvature is circular rather than
parabolic, and the relevant mixing challenge is angular rather than
along a curved ridge. Methods with diverse proposal directions (WMC,
which perturbs along multiple random directions simultaneously) may show
more uniform angular coverage than methods that propose along a single
interpolated line (QMC, AIES).

### Trace plots: ring \\x_1\\

``` r

traces_ring <- do.call(rbind, lapply(diag_methods, function(m) {
  chain <- results_ring[[m]]$chain
  data.frame(
    iteration = seq_len(nrow(chain)),
    x1 = chain[, 1],
    method = m
  )
}))
traces_ring$method <- factor(traces_ring$method, levels = diag_methods)

ggplot(traces_ring, aes(x = iteration, y = x1, color = method)) +
  geom_line(alpha = 0.5, linewidth = 0.3) +
  facet_wrap(~method, ncol = 2) +
  scale_color_manual(values = method_colors) +
  labs(x = "Iteration (post burn-in)", y = expression(x[1]),
       title = "Ring distribution: trace plots") +
  theme_bw() +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))
```

For the ring distribution, the \\x_1\\ trace should oscillate between
approximately \\-5\\ and \\+5\\ as the sampled position moves around the
circle. A trace that stays near \\+5\\ (or any fixed value) indicates
the walker is stuck at one angular position. All methods should achieve
reasonable angular mixing on this 2D problem.

### Autocorrelation: ring \\x_1\\

``` r

max_lag <- 100
acf_ring <- do.call(rbind, lapply(diag_methods, function(m) {
  chain <- results_ring[[m]]$chain[, 1]
  ac <- acf(chain, lag.max = max_lag, plot = FALSE)
  data.frame(lag = as.numeric(ac$lag), acf = as.numeric(ac$acf), method = m)
}))
acf_ring$method <- factor(acf_ring$method, levels = diag_methods)

ggplot(acf_ring, aes(x = lag, y = acf, color = method)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = method_colors) +
  labs(x = "Lag", y = "Autocorrelation",
       title = expression(paste("Ring: ACF of ", x[1]))) +
  theme_bw() +
  theme(legend.position = "bottom")
```

Faster autocorrelation decay means more independent samples per
iteration. Methods whose ACF drops below zero within 10–20 lags are
mixing well; persistent positive autocorrelation beyond 50 lags signals
that the effective sample size is much smaller than the nominal chain
length. On the ring, the \\x_1\\ coordinate oscillates as the walker
traverses the circle, so its autocorrelation structure reflects angular
mixing speed rather than curvature-tracking ability. Some residual
autocorrelation is inherent to ensemble methods because only one
walker’s trajectory is recorded while the ensemble evolves collectively;
the ACF comparison directly explains the ESS differences in the
benchmark tables.

### Posterior recovery: 5D normal

``` r

forest_df <- do.call(rbind, lapply(diag_methods, function(m) {
  r <- results_normal[[m]]
  mcse <- sapply(seq_len(ncol(r$chain)), function(j) MCSE(r$chain[, j]))
  data.frame(
    parameter = paste0("x", 1:5),
    mean = r$post_mean,
    lower = r$post_mean - 1.96 * mcse,
    upper = r$post_mean + 1.96 * mcse,
    method = m
  )
}))
forest_df$method <- factor(forest_df$method, levels = diag_methods)

ggplot(forest_df, aes(x = mean, y = parameter, color = method)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbarh(aes(xmin = lower, xmax = upper),
                 position = position_dodge(width = 0.5), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = method_colors) +
  labs(x = "Posterior mean (95% MCSE interval)", y = NULL,
       title = "5D normal: parameter recovery") +
  theme_bw() +
  theme(legend.position = "bottom")
```

The dashed line shows the true mean of zero. All methods should produce
95% Monte Carlo standard error intervals that contain the true value.
The interval widths reflect each method’s MCSE: narrower intervals
indicate more efficient sampling. On this isotropic target, all methods
should perform comparably, confirming that the QMC methods do not
introduce bias on simple targets.

  

## Practical guidance

### Choosing the right method

The benchmark results point to a clear decision heuristic based on
posterior geometry. The following table summarizes when each method
class is most appropriate.

| Geometry | Recommended | Rationale |
|:---|:---|:---|
| Isotropic / near-Gaussian | WMC or AIES | No curvature to exploit; avoid unnecessary Jacobian penalty |
| Banana / curved ridge | QMC, QMC3, or DQMC | Polynomial proposals track the parabolic/cubic manifold |
| Strong unimodal curvature | DQMC (small b) | Energy-directed proposals target the mode efficiently |
| Ring / angular | WMC or SAMC | Diverse perturbation directions matter more than curvature tracking |
| High-dimensional (N \> 10) | WMC (Gaussian knots if QMC needed) | Jacobian \|w_i\|^N grows exponentially; WMC avoids it entirely |
| Unknown / exploratory | WMC first, then QMC | WMC is robust across geometries; upgrade to QMC if mixing is poor |

Method selection heuristics by posterior geometry {.table}

### Tuning the scale parameter `a`

The `a` parameter is the most important tuning knob for QMC-type methods
(QMC, QMC3, QMCN, DQMC, SQMC). It controls the range of the knot
distribution and hence the spread of Lagrange weights; because the
Jacobian penalty \\\|w_i\|^N\\ grows exponentially with parameter
dimension \\N\\, smaller `a` is needed for higher-dimensional problems.
The benchmark results confirm this scaling: `a = 0.8` works well in 2
dimensions, but `a = 0.3` is needed for 8 dimensions to maintain
reasonable acceptance rates. As a rough guide, start with \\a \approx
1.5/\sqrt{N}\\ and adjust until the acceptance rate falls in the 0.2–0.5
range.

Setting `Gaussian = TRUE` switches from uniform to Gaussian knot
sampling. The Gaussian distribution concentrates mass near \\t = 0\\,
which naturally favors proposals with \\w_i \approx 1\\ (close to the
current position). Militzer [\[4\]](#ref4) recommends Gaussian sampling
for \\N \> 10\\, where uniform sampling wastes too much probability mass
on extreme knots that produce rejected proposals. For \\N \leq 10\\,
uniform sampling is simpler and generally preferable.

For affine methods (MAMC, SAMC), the `a` parameter controls the stretch
range \\\[1/a, a\]\\ and must be greater than 1. The standard AIES
recommendation of \\a = 2\\ works well in most settings. Larger `a`
produces more exploratory moves but with higher rejection rates; the
scaling with dimension is milder than for QMC methods because
\\z^{N-1}\\ is typically less extreme than \\\|w_i\|^N\\ for the same
effective step size.

### Guide points and simplex averaging

The `nGuidePoints` parameter in SQMC, SAMC, and WMC controls how many
walkers contribute to the guide information. More guide points reduce
proposal noise through averaging but require proportionally larger
ensembles (\\N_c \geq 2n_G + 1\\ for SQMC, \\N_c \geq n_G + 1\\ for SAMC
and WMC). The `weightRange` parameter in SQMC and SAMC controls the
randomness of the averaging weights; `weightRange = 1.0` gives uniformly
random weights (full mixing), while `weightRange = 0.0` gives equal
weights (pure centroid).

### Reparameterization

Because all QMC-type proposals involve weighted averages of walker
positions, they implicitly assume that the high-density region is
connected through the affine hull of the ensemble. For parameters with
bounded support (positive parameters, correlations, probabilities),
reparameterization to an unconstrained space is essential: use
log-transform for positive parameters, logit-transform for bounded
parameters, and Cholesky factors for covariance matrices. Without
reparameterization, proposals that average walker positions on the
constrained scale can produce values outside the support, leading to
rejected proposals and wasted computation.

### Computational cost

All methods have similar per-iteration computational cost: one model
evaluation per walker per iteration, giving \\N_c\\ evaluations per
iteration. The model evaluation typically dominates the runtime, so ESS
per model evaluation (equivalently, ESS per iteration divided by
\\N_c\\) is the most relevant efficiency metric for comparing methods.
When the model is expensive, maximizing ESS per evaluation matters more
than maximizing acceptance rate.

### Limitations

Several limitations apply uniformly across the QMC family. First, all
ensemble methods require a population of walkers (\\N_c \geq 3\\ at
minimum, and typically \\2N\\ to \\3N\\ for \\N\\ parameters), so each
iteration costs \\N_c\\ model evaluations; for very expensive models,
this can be prohibitive compared to single-chain methods like NUTS that
achieve high ESS with one evaluation per leapfrog step. Second, the
methods assume continuous, unconstrained parameter spaces; discrete
parameters or hard boundaries require augmentation or
reparameterization. Third, ensemble methods do not exploit gradient
information (except DQMC, which uses energy values but not gradients),
so they cannot match the per-evaluation efficiency of gradient-based
methods on smooth, high-dimensional targets. Their advantage lies in
gradient-free operation, affine invariance, and the ability to track
nonlinear manifold geometry through polynomial interpolation.

  

## Specs reference

| Method |  Min Nc  | Default a | Gaussian |           Extra           |     Jacobian     |
|:-------|:--------:|:---------:|:--------:|:-------------------------:|:----------------:|
| QMC    |    3     |     1     |   opt    |            \-             |    \|w_i\|^N     |
| QMC3   |    4     |     1     |   opt    |            \-             |    \|w_i\|^N     |
| QMCN   | nOrder+1 |     1     |   opt    |          nOrder           |    \|w_i\|^N     |
| DQMC   |    3     |     1     |   opt    |             b             | \|w_i\|^N + asym |
| SQMC   | 2\*nG+1  |     1     |   opt    | nGuidePoints, weightRange |    \|w_i\|^N     |
| MAMC   |    3     |     2     |    \-    |            \-             |     z^(N-1)      |
| SAMC   |   nG+1   |     2     |    \-    | nGuidePoints, weightRange |     z^(N-1)      |
| WMC    |   nG+1   |     1     |    \-    |       nGuidePoints        |       none       |

Summary of QMC ensemble methods {.table}

The `Gaussian` column indicates whether the method supports the
`Gaussian = TRUE` option for Gaussian knot sampling (marked “opt” for
optional). Methods marked “-” use the stretch distribution or symmetric
walk and do not sample knots. The `Jacobian` column shows the acceptance
correction: \\\|w_i\|^N\\ for Lagrange-based methods, \\z^{N-1}\\ for
affine methods, “\|w_i\|^N + asym” for DQMC which adds the asymmetric
forward/reverse correction, and “none” for the symmetric WMC walk.

  

## References

**\[1\]** Goodman, J. and Weare, J. (2010). Ensemble samplers with
affine invariance. *Communications in Applied Mathematics and
Computational Science*, 5(1), 65–80.
[doi:10.2140/camcos.2010.5.65](https://doi.org/10.2140/camcos.2010.5.65)

**\[2\]** Foreman-Mackey, D., Hogg, D.W., Lang, D. and Goodman, J.
(2013). emcee: the MCMC hammer. *Publications of the Astronomical
Society of the Pacific*, 125(925), 306–312.
[doi:10.1086/670067](https://doi.org/10.1086/670067)

**\[3\]** Militzer, B. (2023). Study of Jupiter’s interior with
Quadratic Monte Carlo simulations. *The Astrophysical Journal*, 953,
111.
[doi:10.3847/1538-4357/ace1f1](https://doi.org/10.3847/1538-4357/ace1f1)

**\[4\]** Militzer, B. (2025). Ensemble Monte Carlo calculations with
five novel moves. *Computer Physics Communications*, 307, 109424.
[doi:10.1016/j.cpc.2024.109424](https://doi.org/10.1016/j.cpc.2024.109424)

**\[5\]** ter Braak, C.J.F. (2006). A Markov chain Monte Carlo version
of the genetic algorithm Differential Evolution: easy Bayesian computing
for real parameter spaces. *Statistics and Computing*, 16, 239–249.
[doi:10.1007/s11222-006-8769-1](https://doi.org/10.1007/s11222-006-8769-1)

**\[6\]** Christen, J.A. and Fox, C. (2010). A general purpose sampling
algorithm for continuous distributions (the t-walk). *Bayesian
Analysis*, 5(2), 263–282.
[doi:10.1214/10-BA603](https://doi.org/10.1214/10-BA603)

**\[7\]** Flyvbjerg, H. and Petersen, H.G. (1989). Error estimates on
averages of correlated data. *Journal of Chemical Physics*, 91, 461–466.
[doi:10.1063/1.457480](https://doi.org/10.1063/1.457480)

**\[8\]** Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and
Burkner, P.-C. (2021). Rank-normalization, folding, and localization: an
improved \\\hat{R}\\ for assessing convergence of MCMC. *Bayesian
Analysis*, 16(2), 667–718.
[doi:10.1214/20-BA1221](https://doi.org/10.1214/20-BA1221)

**\[9\]** Geyer, C.J. (1992). Practical Markov chain Monte Carlo.
*Statistical Science*, 7(4), 473–483.
